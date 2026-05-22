# Kiro CLI Troubleshooting — `os error 2`

If you see `error: No such file or directory (os error 2)` when running `kiro-cli chat`, work through these steps in order.

---

## Step 1: Run diagnostics

```bash
kiro-cli diagnostic
```

This prints your version, install method, OS, and PATH. Send the output to John if nothing below fixes it.

## Step 2: Deactivate conda

If you're in a conda environment (you'll see `(base)` or similar in your prompt), conda may be hiding the kiro binary:

```bash
conda deactivate
kiro-cli chat
```

If that fixes it, add this to your `~/.zshrc` to make it permanent:

```bash
export PATH="/opt/homebrew/bin:$PATH"
```

Then restart your terminal.

## Step 3: Check your PATH

```bash
which kiro-cli
echo $PATH
```

`which` should return something like `/opt/homebrew/bin/kiro-cli` or `/usr/local/bin/kiro-cli`. If it returns nothing, kiro-cli isn't installed properly.

## Step 4: Reinstall via Homebrew

```bash
brew reinstall kiro-cli
brew link kiro-cli
```

If you didn't install via Homebrew, download from: https://kiro.dev/downloads/

## Step 5: Try verbose mode

```bash
kiro-cli chat -v
kiro-cli chat -vvv
```

This prints debug output showing exactly which file it's trying to find.

## Step 6: Check the log file

After a failed attempt:

```bash
cat $TMPDIR/kiro-log/kiro-chat.log
```

## Step 7: Check shell and permissions

```bash
echo $SHELL
ls -la $(which kiro-cli)
```

Shell should be `/bin/zsh` or `/bin/bash`. The binary should have execute permissions (`-rwxr-xr-x`).

## Step 8: Clean reinstall

```bash
brew uninstall kiro-cli 2>/dev/null
rm -rf ~/.kiro/sessions
brew install kiro-cli
kiro-cli --version
kiro-cli login --license pro --identity-provider <YOUR_SSO_URL> --region us-east-1
kiro-cli chat
```

---

## Still stuck?

Send John the output of `kiro-cli diagnostic` and the contents of `$TMPDIR/kiro-log/kiro-chat.log`. Happy to screen share — no voice needed, we can troubleshoot over chat.
