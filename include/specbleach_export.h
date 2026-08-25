/*
libspecbleach - A spectral processing library

Copyright 2022 Luciano Dato <lucianodato@gmail.com>

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef SPECBLEACH_EXPORT_H_INCLUDED
#define SPECBLEACH_EXPORT_H_INCLUDED

/**
 * Symbol visibility annotations.
 *
 * Every public library function must be declared with SPECBLEACH_API.
 * Internal symbols must NOT be exported: builds default to hidden
 * visibility so only annotated declarations are visible in the dynamic
 * symbol table.
 *
 * Build system contract:
 * - Building the shared library defines SPECBLEACH_EXPORTS.
 * - Consuming the library as a static archive defines SPECBLEACH_STATIC.
 */
#ifdef SPECBLEACH_STATIC
#define SPECBLEACH_API
#elif defined(_WIN32) || defined(__CYGWIN__)
#ifdef SPECBLEACH_EXPORTS
#define SPECBLEACH_API __declspec(dllexport)
#else
#define SPECBLEACH_API __declspec(dllimport)
#endif
#else
#if defined(__GNUC__) && __GNUC__ >= 4
#define SPECBLEACH_API __attribute__((visibility("default")))
#else
#define SPECBLEACH_API
#endif
#endif

#endif /* SPECBLEACH_EXPORT_H_INCLUDED */
