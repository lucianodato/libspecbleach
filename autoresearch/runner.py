import os
import sys
import subprocess
import json
import time
from dotenv import load_dotenv
from pydantic import BaseModel, Field
from google import genai
from google.genai import types

MODEL_TEMPERATURE = 0.4
MAX_ITERATIONS = 10
BEST_SCORE = -999.0

# Priority list using your specified models
MODEL_FALLBACK_LIST = [
    "gemini-3.6-flash",
    "gemini-3.5-flash",
    "gemini-3.5-flash-lite"
]

# Anchor working directory to repository root
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR) if os.path.basename(SCRIPT_DIR) == "autoresearch" else SCRIPT_DIR
os.chdir(REPO_ROOT)

# --- LOAD API KEY FROM .env FILE ---
dotenv_path_autoresearch = os.path.join(SCRIPT_DIR, ".env")
dotenv_path_root = os.path.join(REPO_ROOT, ".env")

if os.path.exists(dotenv_path_autoresearch):
    load_dotenv(dotenv_path_autoresearch)
elif os.path.exists(dotenv_path_root):
    load_dotenv(dotenv_path_root)
else:
    load_dotenv()

if not os.environ.get("GEMINI_API_KEY"):
    print("[CRITICAL ERROR] GEMINI_API_KEY not found in environment or .env file!")
    print("Please create a .env file containing: GEMINI_API_KEY=your_key_here")
    sys.exit(1)

client = genai.Client()

class CEditProposal(BaseModel):
    filepath: str = Field(description="Relative path to C source file, e.g., 'src/denoiser_core.c'")
    new_code: str = Field(description="The full updated content of the C file")
    hypothesis: str = Field(description="A 1-sentence explanation of the DSP change")

def get_allowed_c_code() -> str:
    """Reads C source and header files from root src/ and include/ directories."""
    code_str = ""
    for folder in ["src", "include"]:
        if not os.path.exists(folder):
            continue
        for root, _, files in os.walk(folder):
            for f in files:
                if f.endswith(".c") or f.endswith(".h"):
                    path = os.path.join(root, f)
                    with open(path, "r", encoding="utf-8") as file:
                        code_str += f"\n--- FILE: {path} ---\n" + file.read()
    return code_str

def generate_patch_with_fallback(prompt: str) -> CEditProposal:
    """Tries generating code across your specified models, failing over automatically on quota limits."""
    last_error = None
    for model_name in MODEL_FALLBACK_LIST:
        try:
            print(f"Querying model [{model_name}]...")
            
            config = types.GenerateContentConfig(
                response_mime_type="application/json",
                response_schema=CEditProposal,
                temperature=MODEL_TEMPERATURE,
                thinking_config=types.ThinkingConfig(thinking_budget=1024)
            )

            response = client.models.generate_content(
                model=model_name,
                contents=prompt,
                config=config,
            )
            
            # 1. Try SDK's structured Pydantic parser
            if getattr(response, "parsed", None):
                return response.parsed
                
            # 2. Manual JSON fallback if parsed is None
            if getattr(response, "text", None):
                text = response.text
                clean_json = text[text.find("{"):text.rfind("}") + 1]
                data = json.loads(clean_json)
                return CEditProposal(**data)

        except Exception as e:
            last_error = e
            print(f"[WARNING] Model [{model_name}] failed: {e}")
            print("-> Switching to next model in fallback list...\n")
            time.sleep(2)
            continue

    raise RuntimeError(f"All models in MODEL_FALLBACK_LIST failed! Last error: {last_error}")

def run_evaluator() -> float:
    """Executes evaluator.py in FAST mode (6 cases, 1 preset)."""
    eval_script = os.path.join(SCRIPT_DIR, "evaluator.py")
    cmd = [sys.executable, eval_script, "--fast", "--max-cases", "6"]
    
    res = subprocess.run(cmd, capture_output=True, text=True)
    for line in res.stdout.split("\n"):
        if "METRIC_SCORE:" in line:
            try:
                return float(line.split(":")[1].strip())
            except ValueError:
                pass
    return -999.0

def main():
    global BEST_SCORE
    
    print(f"Root Repository Directory: {os.getcwd()}")
    print("Evaluating baseline C codebase...")
    BEST_SCORE = run_evaluator()
    print(f"Initial Baseline Score: {BEST_SCORE:.4f}\n")

    if BEST_SCORE == -999.0:
        print("[CRITICAL ERROR] Baseline build/evaluation failed. Check CMake build manually!")
        sys.exit(1)

    for i in range(MAX_ITERATIONS):
        print(f"=================== ITERATION {i+1} ===================")
        c_code = get_allowed_c_code()
        
        prompt = f"""You are an expert Audio DSP C Developer optimizing libspecbleach.

Read program.md carefully. Your task is to propose STRUCTURAL algorithmic changes to the C processing logic in `src/`.

Current Best Composite Score: {BEST_SCORE:.4f}

CRITICAL RULES & CONSTRAINTS:
1. DO NOT use heavy recursive frame-to-frame smoothing, Decision-Directed SNR tracking with high alpha, or spectral delay loops that cause ECHOES and temporal ghosting!
2. DO PREFER 2D Image Processing techniques on the spectrogram (e.g., Bilateral filtering, Guided filtering, Anisotropic diffusion, 2D Morphological filtering) or instantaneous psychoacoustic masking.
3. DO NOT simply tweak constant numbers or threshold multipliers.
4. Ensure your C code compiles cleanly via CMake and manages memory safely without leaks.

C CODEBASE CONTEXT:
{c_code[:100000]}
"""

        filepath = None
        try:
            patch: CEditProposal = generate_patch_with_fallback(prompt)
            filepath = patch.filepath

            # --- GUARDRAIL: Prevent LLM from editing anything outside src/ or include/ ---
            clean_path = os.path.normpath(filepath)
            if not (clean_path.startswith("src") or clean_path.startswith("include")):
                print(f"[REJECTED] Agent attempted to modify a forbidden path: {filepath}")
                continue

            os.makedirs(os.path.dirname(filepath), exist_ok=True)
            
            print(f"Hypothesis: {patch.hypothesis}")
            print(f"Applying edit to {filepath}...")

            with open(filepath, "w", encoding="utf-8") as f:
                f.write(patch.new_code)

            new_score = run_evaluator()
            print(f"New Score: {new_score:.4f} vs Best Score: {BEST_SCORE:.4f}")

            if new_score > BEST_SCORE:
                BEST_SCORE = new_score
                print(">>> IMPROVEMENT ACCEPTED! Committing to git...")
                subprocess.run(["git", "add", filepath])
                subprocess.run(["git", "commit", "-m", f"autoresearch: {patch.hypothesis} (score: {new_score:.4f})"])
            else:
                print("<<< REJECTED / CRASHED. Rolling back target file...")
                # Revert ONLY the specific C file edited during this step
                subprocess.run(["git", "checkout", "HEAD", "--", filepath])

        except Exception as e:
            print(f"Error in iteration {i+1}: {e}")
            if filepath and os.path.exists(filepath):
                # Revert specific file if known
                subprocess.run(["git", "checkout", "HEAD", "--", filepath])
            else:
                # Revert ONLY src/ and include/ directories on general crash—never autoresearch/
                subprocess.run(["git", "checkout", "HEAD", "--", "src/", "include/"])
            
        time.sleep(1)

if __name__ == "__main__":
    main()
