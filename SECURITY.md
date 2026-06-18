# Security Policy

## Supported Versions

| Version | Supported |
|---------|-----------|
| 0.2.x   | ✅        |
| < 0.2   | ❌        |

## Reporting a Vulnerability

If you discover a security vulnerability, please report it privately:

**Email:** john@nodes.bio

Do **not** open a public GitHub issue for security vulnerabilities.

We will acknowledge receipt within 48 hours and provide a timeline for a fix.

## Scope

Rosetta executes R code via rpy2. Security concerns related to:
- Arbitrary code execution through crafted design formulas
- Path traversal in file-based inputs
- R package supply chain issues

are all in scope. Please report anything that feels off.
