from fractions import Fraction
import math
from functools import reduce


def poly_eval(poly, x):
	"""Evaluate polynomial with coefficients 'poly' (in descending order) at x."""
	result = Fraction(0)
	for coeff in poly:
		result = result * x + coeff
	return result


def synthetic_division(poly, r):
	"""
	Perform synthetic division of a polynomial by (x - r).
	poly: list of coefficients (as Fractions) in descending order.
	r: the root (Fraction) such that (x - r) is a factor.
	Returns (quotient, remainder).
	"""
	n = len(poly)
	quotient = [poly[0]]
	for i in range(1, n):
		quotient.append(poly[i] + r * quotient[i - 1])
	remainder = quotient.pop()  # last element is the remainder
	return quotient, remainder


def divisors(n):
	"""Return a set of all (positive and negative) divisors of integer n."""
	n = abs(n)
	divs = set()
	for i in range(1, int(math.sqrt(n)) + 1):
		if n % i == 0:
			divs.add(i)
			divs.add(n // i)
	# Include negatives as well
	all_divs = set()
	for d in divs:
		all_divs.add(d)
		all_divs.add(-d)
	return all_divs


def possible_rational_roots(poly):
	"""
	Generate all possible rational roots for a polynomial with integer coefficients.
	Based on factors of the constant term and the leading coefficient.
	"""
	lc = poly[0]
	ct = poly[-1]
	p_divs = divisors(ct)
	q_divs = divisors(lc)
	candidates = set()
	for p in p_divs:
		for q in q_divs:
			candidates.add(Fraction(p, q))
	return candidates


def factor_quadratic(poly):
	"""
	Factor a quadratic polynomial poly = [a, b, c] with a, b, c integers.
	If the roots are rational, return the factors in the form (A, B) representing A*x+B.
	Otherwise, return the quadratic as irreducible.
	"""
	a, b, c = poly
	D = b * b - 4 * a * c
	# Check for perfect square discriminant
	sqrtD = int(math.isqrt(D))
	if sqrtD * sqrtD != D:
		return [poly]  # irreducible quadratic over the rationals
	# Compute roots using the quadratic formula
	root1 = Fraction(-b + sqrtD, 2 * a)
	root2 = Fraction(-b - sqrtD, 2 * a)

	def factor_from_root(root):
		# Write the factor as (denom*x - num) so that coefficients are integers.
		return (root.denominator, -root.numerator)

	return [factor_from_root(root1), factor_from_root(root2)]


def factor_polynomial(poly):
	"""
	Factor a polynomial with integer coefficients using a Kronecker-like method.
	poly: list of coefficients (in descending order).
	Returns a tuple (constant_factor, factors) where 'factors' is a list.
	Each factor is either a tuple (A, B) representing (A*x + B) or
	a list of coefficients for an irreducible polynomial factor.
	"""

	# Remove common content (gcd of coefficients) and make the leading coefficient positive.
	def gcd_list(nums):
		return reduce(math.gcd, nums)

	cont = gcd_list([abs(c) for c in poly])
	constant_factor = 1
	if poly[0] < 0:
		poly = [-c for c in poly]
		constant_factor *= -1
	if cont != 0:
		poly = [c // cont for c in poly]
		constant_factor *= cont

	def factor_rec(poly):
		n = len(poly) - 1  # degree of the polynomial
		if n == 1:
			# Linear polynomial: a*x + b
			a, b = poly
			return [(a, b)]
		elif n == 2:
			# Try to factor quadratic
			quad_factors = factor_quadratic(poly)
			if len(quad_factors) == 1:  # irreducible
				return [poly]
			else:
				return quad_factors
		else:
			# Look for a rational root using the possible candidates.
			candidates = possible_rational_roots(poly)
			for r in candidates:
				if poly_eval([Fraction(c) for c in poly], r) == 0:
					# Found a rational root r.
					# The corresponding linear factor is (x - r). To have integer coefficients,
					# write it as (q*x - p) where r = p/q in lowest terms.
					factor_factor = (r.denominator, -r.numerator)
					# Perform synthetic division to obtain the quotient.
					q, rem = synthetic_division([Fraction(c) for c in poly], r)
					# Clear denominators in quotient.
					lcm = 1
					for frac in q:
						lcm = lcm * frac.denominator // math.gcd(lcm, frac.denominator)
					q_int = [int(frac * lcm) for frac in q]
					return [factor_factor] + factor_rec(q_int)
			# If no factor is found, return the polynomial as irreducible.
			return [poly]

	factors = factor_rec(poly)
	return constant_factor, factors


# Define the polynomial coefficients for -4x^3 - 6x^2 + 58x + 120.
poly = [-4, -6, 58, 120]
cf, factors = factor_polynomial(poly)


# Format the factors nicely.
def format_factor(f):
	if isinstance(f, tuple):
		A, B = f
		sign = "+" if B >= 0 else "-"
		B_abs = abs(B)
		# If A is 1, we write just x; otherwise, include A.
		A_str = "" if A == 1 else str(A)
		return f"({A_str}x {sign} {B_abs})"
	else:
		# Irreducible polynomial factor given by its coefficients.
		return str(f)


factor_strs = [format_factor(f) for f in factors]
result = f"{cf} " + " * ".join(factor_strs)
print("The factorization is:")
print(result)
