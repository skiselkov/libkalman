/*
 * CDDL HEADER START
 *
 * The contents of this file are subject to the terms of the
 * Common Development and Distribution License, Version 1.0 only
 * (the "License").  You may not use this file except in compliance
 * with the License.
 *
 * You can obtain a copy of the license in the file COPYING
 * or http://www.opensource.org/licenses/CDDL-1.0.
 * See the License for the specific language governing permissions
 * and limitations under the License.
 *
 * When distributing Covered Code, include this CDDL HEADER in each
 * file and include the License file COPYING.
 * If applicable, add the following below this CDDL HEADER, with the
 * fields enclosed by brackets "[]" replaced with your own identifying
 * information: Portions Copyright [yyyy] [name of copyright owner]
 *
 * CDDL HEADER END
 */
/*
 * Copyright 2020 Saso Kiselkov. All rights reserved.
 */

#ifndef	_KALMAN_ASSERT_H_
#define	_KALMAN_ASSERT_H_

#include <assert.h>
#include <stdlib.h>

#ifdef	KALMAN_USE_ACFUTILS
#include <acfutils/safe_alloc.h>
#endif

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef	KALMAN_USE_ACFUTILS

#define	KAL_ASSERT(x)		ASSERT(x)
#define	KAL_ASSERT3U(x, y, z)	ASSERT3U(x, y, z)
#define	KAL_ASSERT3F(x, y, z)	ASSERT3F(x, y, z)
#define	KAL_VERIFY(x)		VERIFY(x)

#else	/* !defined(KALMAN_USE_ACFUTILS) */

static inline void *
safe_calloc(size_t nmemb, size_t size)
{
	void *p = calloc(nmemb, size);
	if (nmemb > 0 && size > 0)
		assert(p != NULL);
	return (p);
}
#define	KAL_ASSERT(x)		assert(x)
#define	KAL_ASSERT3U(x, y, z)	assert(x y z)
#define	KAL_ASSERT3F(x, y, z)	assert(x y z)
#define	KAL_VERIFY(x)		assert(x)

#if	__STDC_VERSION__ >= 201112L
#define	CTASSERT(x)	_Static_assert((x), #x)
#else	/* __STDC_VERSION__ < 201112L */
#if	defined(__GNUC__) || defined(__clang__)
#define	CTASSERT(x)		_CTASSERT(x, __LINE__)
#define	_CTASSERT(x, y)		__CTASSERT(x, y)
#define	__CTASSERT(x, y)	\
	typedef char __compile_time_assertion__ ## y [(x) ? 1 : -1] \
	    __attribute__((unused))
#else	/* !defined(__GNUC__) && !defined(__clang__) */
#define	CTASSERT(x)
#endif	/* !defined(__GNUC__) && !defined(__clang__) */
#endif	/* __STDC_VERSION__ < 201112L */

#endif	/* !defined(KALMAN_USE_ACFUTILS) */

#ifdef	__cplusplus
}
#endif

#endif	/* _KALMAN_ASSERT_H_ */
