# Lessons

- When a user asks to re-check code after changes, do not rely on prior findings. Re-inspect the current repository state first, preferably by reviewing the latest commit/diff before reporting.
- For interactive CLI defaults, do not assume CSV-on-stdout is acceptable for exploratory commands. If the user is expected to run a long scientific workflow directly, provide visible progress and a readable default report, while keeping machine-readable modes explicit.
