# Lessons

- When a user asks to re-check code after changes, do not rely on prior findings. Re-inspect the current repository state first, preferably by reviewing the latest commit/diff before reporting.
- For interactive CLI defaults, do not assume CSV-on-stdout is acceptable for exploratory commands. If the user is expected to run a long scientific workflow directly, provide visible progress and a readable default report, while keeping machine-readable modes explicit.
- When a progress complaint includes a traceback, inspect the interrupted stack for the actual bottleneck before adding cosmetic logging. Prefer changing the slow call path when a cheaper equivalent exists.
- When a console script can run but an import fails, check the script shebang and `sys.executable` before assuming the package is absent. CLI availability on `PATH` may come from a different virtualenv than the command being debugged.
- When explaining CTA panel gates, distinguish the original adaptive HPA restriction filter from later panel-specific safety vetoes. Report the exact source columns and installed CLI import path/version before concluding a selector default is behaving as intended.
