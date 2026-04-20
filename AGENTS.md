## Golden Rules

1. **Never commit to `main`.** Always `git checkout -b <feature-branch>` before editing. Land via PR.
2. **Every PR bumps the version.** Even doc-only PRs — at minimum a patch bump. The version lives in `pyensembl/version.py` and must be bumped in the PR itself.
3. **"Done" means merged AND deployed to PyPI** — never stop at merge. After a PR merges, run `./deploy.sh` from a clean `main`. Skipping deploy = task not done.
4. **File problems as issues, don't silently work around them.** If you hit a bug here or in a sibling openvax repo, open a GitHub issue on the correct repo and link it from the PR.
5. **After a PR ships, look for the next block of work.** Read open issues across the relevant openvax repos, group by dependency + urgency. Prefer *foundational* changes that unblock multiple downstream improvements; otherwise chain the smallest independent improvements.

---

## Before Completing Any Task

Before considering any code change complete, you MUST:

1. **Run `./lint.sh`** — runs `ruff check pyensembl/`.
2. **Run `./test.sh`** — runs the pytest suite.

Do not tell the user you are "done" or that changes are "complete" until both pass. `./lint-and-test.sh` runs them together.

## Scripts

- `./lint.sh` — runs `ruff check pyensembl/`. **Always use this for linting.**
- `./test.sh` — runs pytest with coverage (must pass).
- `./lint-and-test.sh` — convenience wrapper that runs lint then tests.
- `./deploy.sh` — deploys to PyPI. Gates on `lint.sh` + `test.sh`, builds sdist+wheel, uploads via twine, tags the commit with `pyensembl/version.py`, pushes tags. **Bump `pyensembl/version.py` before running** — deploy.sh does not bump the version for you.
- `./develop.sh` — installs package in development mode.

## Code Style

- Use ruff for linting (there is no `format.sh`; formatting is not currently enforced in CI).
- Configuration is in `pyproject.toml`.
- Python support: 3.9+.

---

## Workflow Orchestration

### 1. Upfront Planning
- For ANY non-trivial task (3+ steps or architectural decisions): write a detailed spec before touching code
- If something goes sideways, STOP and re-plan immediately — don't keep pushing
- Use planning/verification steps, not just building
- Write detailed specs upfront to reduce ambiguity

### 2. Self-Improvement Loop
- After ANY correction from the user: capture the pattern (in Claude Code memory or `tasks/lessons.md`)
- Write rules for yourself that prevent the same mistake
- Ruthlessly iterate on these lessons until mistake rate drops
- Review lessons at session start for relevant project

### 3. Verification Before Done
- Never mark a task complete without proving it works
- Diff behavior between the latest code and your changes when relevant
- Ask yourself: "Would a staff engineer approve this?"
- Run tests, check logs, demonstrate correctness

### 4. Demand Elegance (Balanced)
- For non-trivial changes: pause and ask "is there a more elegant way?"
- If a fix feels hacky: "Knowing everything I know now, implement the elegant solution"
- Skip this for simple, obvious fixes — don't over-engineer
- Challenge your own work before presenting it

### 5. Autonomous Bug Fixing
- When given a bug report: just fix it. Don't ask for hand-holding
- Point at logs, errors, failing tests — then resolve them
- Zero context switching required from the user
- Fix failing unit tests without being told how

---

## Core Principles

- **Simplicity First**: Make every change as simple as possible. Impact minimal code.
- **No Laziness**: Find root causes. No temporary fixes. Senior developer standards.
- **Minimal Impact**: Changes should only touch what's necessary. Avoid introducing bugs.
- **No tautological tests**: Don't write tests that reassert the contents of declarative config (e.g. a `pyproject.toml` dependency list against a hardcoded copy). They verify nothing and break on every legitimate bump.

## Scientific Domain Knowledge
- **Read the literature**: if some code involves scientific or biological concepts, feel free to search for review papers and read those before changing code that expresses scientific concepts.
- **Flag inconsistencies**: if code expresses a scientific model that's at odds with your understanding, note that inconsistency and ask for clarification.
