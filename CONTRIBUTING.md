# The RangeShifter platform - An eco-evolutionary modelling framework

## How to contribute

Thank you for your interest in contributing to the RangeShifter platform. 
In this document we will give you guidance on how to contribute to the RangeShifter project regarding issues, bug fixing and adding new features. In this guidance we distinguish between contributing to the RangeShifter core code and the R interface (RangeShiftR package).

## Repo structure
![](https://raw.githubusercontent.com/RangeShifter/RangeShiftR-package-dev/development-guidelines/RangeShiftR/man/figures/RS_repos.png?token=GHSAT0AAAAAACKSYLFW67K4UCXY6TPS4WISZK4RXHQ)

RangeShifter is distributed with three user interfaces, each living in their own repo:

- the RangeShifter GUI (clickable Windows interface)*
- RangeShifter Batch Mode (command line interface)
- the RangeShiftR package (R interface)

All three share the same source code for the core simulation (i.e., the actual model), which lives in its own repo (RScore). Each of the interfaces keeps a copy of this core code in a subfolder called RScore, kept in sync with the RScore repo via a git subtree (see Git subtree usage section). 

⚠️ If you wish to propose a change to the core code of the simulation, please do so *in the RScore repo*, rather than in the RScore folder in one of the interfaces. (see below)

*The RangeShifter GUI is currently being rewritten, and is not open source yet.

### RangeShifter core code

As mentioned above, the RangeShifter core code is located and maintained in a dedicated repository RScore, which is currently only maintained by the RangeShifter development team. Any changes regarding the RangeShifter core code should be done in this repository and afterwards synced with all interfaces using the git subtree feature (see [Git subtree usage section](https://github.com/RangeShifter/RangeShiftR-package-dev/blob/development-guidelines/CONTRIBUTING.md#git-subtree-usage)). 

#### Bugs

To report a bug, please [open an issue](https://github.com/RangeShifter/RangeShiftR-package-dev/issues/new), using the Bug Report template. 
Please do check if a related issue has already open on one of the other interfaces ([here](https://github.com/RangeShifter/RangeShifter_batch/issues) for the batch interface).
To propose a bug fix (thank you!!), please create and work on your own branch or fork, from either `main` or `develop` (preferred), and open a pull request when your fix is ready to be merged into the original branch.

**For RangeShifter-batch only, (for now?):** as a prerequisite for merging, please ensure that your version passes status check (that is, RangeShifter can still build and run as intended).
This can be seen in the Actions panel for every commit and at the bottom of the pull request.

Maintainers will review the pull request, possibly request changes, and eventually integrate the bug fix into RScore, and update the subtrees to bring the fix to all interfaces.

#### New features

Do you have an idea of a new feature in the RangeShifter platform that should be integrated and is of use for other RangeShifter users? 
Please get in touch with the RangeShifter development team (rangeshiftr@uni-potsdam.de (*or other mail?*) to discuss a collaboration.

Alternatively*, proceed as with the bug fix above: create your own branch or fork _from `develop`_ and work from there, and submit a pull request when your new features are ready to join the core code. 
We recommend that you update your branch regularly to new changes on `develop` (using `git merge develop`) to reduce the risk of merge conflicts or your version getting out-of-touch in the late stages of development.
We also recommend that you work in small commits, as this makes the code easier to debug, and makes it easier for maintainers to understand your contributions when reviewing a pull request.

*Do we welcome independent contributions?

### RangeShiftR package

#### Issues

Issues should be used for reporting technical problems with the R package or suggest improvements e.g. in the documentation. Questions regarding the application and settings of the RangeShiftR should be asked in the [discussion forum](https://github.com/RangeShifter/RangeshiftR-tutorials/discussions).

##### Create a new issue

If you encounter a technical problem in the RangeShiftR package or would like to suggest an improvement to the package, please search if [a related issue already exists](https://github.com/RangeShifter/RangeShiftR-package-dev/issues). If you can't find a related issue, you can open a new issue using a relevant [issue form](https://github.com/RangeShifter/RangeShiftR-package-dev/issues/new/choose).

##### Solve an issue

As a general rule, the RangeShiftR development team will take care of opened issues. However, feel free to scan through [existing issues](https://github.com/RangeShifter/RangeShiftR-package-dev/issues) and leave comments. 

#### Fixing a bug

If you found a bug in the RangeShiftR package files, which isn't already reported in an [existing issue](https://github.com/RangeShifter/RangeShiftR-package-dev/issues), please report the bug using the dedicated [issue form](https://github.com/RangeShifter/RangeShiftR-package-dev/issues/new/choose). If you were able to solve the bug already, you can suggest the fix in the issue discussion and create a pull request. The pull request will then be reviewed by the RangeShiftR development team and eventually be merged.

#### Adding a new feature

Any suggestions how to improve the R package? Or do you think, we missed a feature in the R package interface? Please check if someone else [already suggested the feature](https://github.com/RangeShifter/RangeShiftR-package-dev/issues) and use the [dedicated form](https://github.com/RangeShifter/RangeShiftR-package-dev/issues/new/choose) to suggest a new feature. You can also contact the RangeShiftR development team directly (rangeshiftr@uni-potsdam.de).

*TODO - move the two sections below into their own document?*

## Git usage

## Git subtree usage

In order to ensure that the same version of RangeShifter's core code is used by all three interfaces (RangeShiftR, RangeShifter-batch and the GUI), each interface repo keeps a copy of RScore as a git subtree. In this section we describe how to interact with RScore via these git subtrees.

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

We haven't yet found a way to push new changes made in a RScore subfolder back into the RScore repo. This is why we ask that contributions are made directly inside the RScore repo.

If you know how to do this, please drop us a line!

Alternatively, if you have already made changes on the subfolder, you could copy its contents into a new branch in RScore, then open a pull request.

### Switching the subfolder to a new branch

There is unfortunately to do so. To track a different branch of RScore, one must delete the RScore subfolder (via git) and import the subtree again:

```bash
git rm src/RScore -r
git commit -m "switching subtree branch"
git subtree add --prefix src/RScore RScore <the new branch>
```


