#!/usr/bin/env bash
# Install the commit-msg hook that enforces the authorship and identifier
# properties in src/check_authorship.py before a message can be recorded.
# Hooks are not version controlled, so run once per clone:  make hooks
set -euo pipefail
root="$(git rev-parse --show-toplevel)"
hook="$root/.git/hooks/commit-msg"
cat > "$hook" <<'HOOK'
#!/usr/bin/env bash
# Refuse a commit message that violates the repository authorship and
# identifier properties. See src/check_authorship.py.
set -euo pipefail
root="$(git rev-parse --show-toplevel)"
python3 -I - "$root" "$1" <<'PY'
import importlib.util, sys
root, msg_file = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location("g", root + "/src/check_authorship.py")
g = importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
hits = g.check_message(open(msg_file, encoding="utf-8").read())
if hits:
    print("commit-msg hook: message refused by the authorship gate")
    for h in hits:
        print("  " + h)
    sys.exit(1)
PY
HOOK
chmod +x "$hook"
echo "installed $hook"
