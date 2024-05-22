# RangeShifter core code 

This repo contains the core simulation code for RangeShifter v2.0 and is not meant to be compiled or run on its own.

<img src="https://github.com/RangeShifter/RScore/blob/main/RScore_logo.png" align="right" height = 200/>

If you are only interested in using RangeShifter, you can ignore this and head to the repo of one of the interfaces:

- [WIP] RangeShifter GUI

- [RangeShiftR](https://github.com/RangeShifter/RangeShiftR-pkg)

- [RangeShifter-batch](https://github.com/RangeShifter/RangeShifter_batch)

## Usage: git subtree

In order to ensure that the same version of RangeShifter's core code is used by all three interfaces (RangeShiftR, RangeShifter-batch and the GUI), each interface repo keeps a copy of RScore as a git subtree. In this section, we describe how to use the git subtrees to update the subfolder and copy this repo anew.

First, in a local clone of one of the interface repos, add a remote named `RScore` pointing to the RScore repo. This will be convenient as a shortcut for git subtree commands.

```bash
git remote add RScore https://github.com/RangeShifter/RScore.git
```

### Pulling new changes

To update the RScore subfolder with new changes made to the RScore repo, one can use the `git subtree pull` command:

```bash
git subtree pull --prefix <path_to_RScore_subfolder> RScore <branch>
```

Note the path must match the location of the RScore subfolder, and the branch must match the one the subtree was originally added from (by default, this should be `main`).

e.g. for RangeShifter-batch, use:

```bash
git subtree pull --prefix src/RScore RScore main
```

while for RangeShiftR, use:

```bash
git subtree pull --prefix RangeShiftR/src/RScore RScore main
```

### Pushing new changes to RScore

```bash
git subtree push --prefix <path_to_RScore_subfolder> RScore <branch>
```

e.g., from RangeShifter-batch's `main` to RScore's `main`:

```bash
git subtree push --prefix src/RScore RScore main
```

### Switching the subfolder to a new branch

There is unfortunately to do so. To track a different branch of RScore, one must delete the RScore subfolder (via git) and import the subtree again:

```bash
git rm src/RScore -r
git commit -m "switching subtree branch"
git subtree add --prefix src/RScore RScore <the new branch>
```

## Contributing to RangeShifter core code

See [CONTRIBUTING](https://github.com/RangeShifter/RScore/blob/main/CONTRIBUTING.md).

## Maintainers

- [@JetteReeg](https://github.com/JetteReeg)
- [@TheoPannetier](https://github.com/TheoPannetier)
