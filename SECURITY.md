# Security policy

Report suspected vulnerabilities privately through GitHub Security Advisories
for `e-south/dense-arrays`. Do not open a public issue containing credentials,
private sequences, or exploit details.

Dense Arrays treats realized-array JSON as untrusted input. Parsers reject
unknown fields and wrong primitive types, and renderers must preserve the
authority and reveal semantics encoded in the playback plan.
