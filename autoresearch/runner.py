import os
import sys
import subprocess
import json
import time
from dotenv import load_dotenv
from pydantic import BaseModel, Field
from google import genai
from google.genai import types

# Anchor working directory to repository root
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR) if os.path.basename(SCRIPT_DIR) == "autoresearch" else SCRIPT_DIR
os.chdir(REPO_ROOT)

# --- LOAD API KEY FROM .env FILE ---
# Check both autoresearch/.env and repo root .env
dotenv_path_autoresearch = os.path.join(SCRIPT_DIR, ".env")
dotenv_path_root = os.path.join(REPO_ROOT, ".env")

if os.path.exists(dotenv_path_autoresearch):
    load_dotenv(dotenv_path_autoresearch)
elif os.path.exists(dotenv_path_root):
    load_dotenv(dotenv_path_root)
else:
    load_dotenv()  # Fallback to standard environment lookup

# Verify API key is loaded
if not os.environ.get("GEMINI_API_KEY"):
    print("[CRITICAL ERROR] GEMINI_API_KEY not found in environment or .env file!")
    print("Please create a .env file containing: GEMINI_API_KEY=your_key_here")
    sys.exit(1)

client = genai.Client()
MODEL_NAME = "gemini-3.6-flash"
MAX_ITERATIONS = 10
BEST_SCORE = -999.0

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

def run_evaluator() -> float:
    """Executes evaluator.py via current python executable."""
    eval_script = os.path.join(SCRIPT_DIR, "evaluator.py")
    res = subprocess.run([sys.executable, eval_script], capture_output=True, text=True)
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
        
        prompt = f"""You are an expert Audio DSP C Developer inventing architectural improvements for libspecbleach.

                    Read program.md carefully. Your task is to propose STRUCTURAL algorithmic changes to the C processing logic in `src/`.

                    Current Best Composite Score: {BEST_SCORE:.4f}

                    CRITICAL RULES:
                    1. DO NOT simply tweak constant numbers, threshold multipliers, or existing parameter scales.
                    2. DO propose structural logic changes (e.g., Multi-Resolution STFT blending, phase-locking, dynamic multi-window masking, or new time-frequency filter topologies).
                    3. If your proposal only changes scalar numbers or magic constants, IT WILL BE REJECTED.

                    Formulate a clear DSP hypothesis, write the complete modified C file, and ensure it builds without warnings or pointer crashes.

                    C CODEBASE CONTEXT:
                    {c_code[:100000]}
                    """


        try:
            response = client.models.generate_content(
                model=MODEL_NAME,
                contents=prompt,
                config=types.GenerateContentConfig(
                    response_mime_type="application/json",
                    response_schema=CEditProposal,
                    temperature=0.3,
                    thinking_config=types.ThinkingConfig(thinking_budget=1024)
                ),
            )
            
            patch: CEditProposal = response.parsed
            filepath = patch.filepath
            
            # Ensure target directory exists
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
                print("<<< REJECTED / CRASHED. Rolling back with git checkout...")
                subprocess.run(["git", "checkout", "HEAD", "--", filepath])

        except Exception as e:
            print(f"Error in iteration {i+1}: {e}")
            subprocess.run(["git", "reset", "--hard", "HEAD"])
            
        time.sleep(1)

if __name__ == "__main__":
    main()
