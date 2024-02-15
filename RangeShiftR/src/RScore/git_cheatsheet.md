# Git Cheatsheet

Quick reference on Git usage for RangeShifter contributors

#### Creating a local copy of this repo

```bash
git clone https://github.com/RangeShifter/RScore.git
```

#### Enquire about the current state of changes

```
git status
```
This will display whether the local branch is up-to-date with its GitHub counterpart, what files have been changed and if they are staged for commit or not.

#### Updating the active branch with latest changes (pulling)

```bash
git pull
```

#### Updating the active branch with another branch (merging)

```bash
git merge <branch to get new changes from>
```

Merging may trigger a merge conflict is the same lines have been changed on both branches. See Solving Merge Conflict below.

#### Staging changes to be committed
Changes to local files must first be added to the commit queue ("staged") before being committed.

```bash
git add <name of file> # single file
git add . # entire active folder
```

#### Creating a commit

```bash
git commit -m "commit message"
```

A good commit message is concise, but descriptive.

#### Uploading local commits to GitHub (pushing)

```bash
git push
```

#### Switching to an existing branch

```bash
git checkout <branch name>
```
Switching does not update the new active branch with GitHub automatically, so make sure to pull after switching!

#### Creating a new branch from active branch

```bash
git branch <name of new branch>
```
You will also need to set the corresponding branch on `origin` (GitHub) before you can push:

```bash
git push --set-upstream origin <name of branch>
```

New branches can also be created on GitHub (drop-down button in the top-left corner of the main page). 
New branches on GitHub are brought to the local copy with the next pull.

#### Deleting a branch locally

```bash
git branch -d <name of branch>
```

#### Instruct Git to not track some files

Open `.gitignore` and add the path to the files or folders to exclude from the git history.
Check with `git status` that the files are indeed not tracked by git anymore.

#### Solving a merge conflict
Merge conflicts can arise when multiple contributors simulatneously change the same line of code.
In such cases, git is incapable of deciding which version should be kept and asks for human input.
Git tells you which files are affected by the conflict.
Open each file and resolve **each** section that looks like this:

```
After opening the box, we can affirm that
<<<<<<<<<<< HEAD               # delete this line
Shroedinger's cat is alive.    # delete either this...
================               # delete this line
Shroedinger's cat is dead.     # ... or this (or keep both, or none, or a different solution)
>>>>>>>>>>> SHA                # delete this line
What an insightful result!
```

Ctrl+F "HEAD" is really helpful for finding the conflicts in large files.
When you are done, create a commit stating e.g. "solved merge conflict" and push.

## Git subtrees

See [Git subtrees](https://github.com/RangeShifter/RScore/tree/development-guidelines#usage-git-subtree) section in README.
