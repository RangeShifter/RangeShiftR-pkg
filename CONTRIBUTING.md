# The RangeShifter platform - An eco-evolutionary modelling framework

## How to contribute

Thank you for your interest in contributing to the RangeShifter platform. 
In this document we will give you guidance on how to contribute to the RangeShifter project regarding issues, bug fixing and adding new features. In this guidance we distinguish between contributing to the RangeShifter core code and the R interface (RangeShiftR package).

## Repo structure
![Rangeshifter repo structure](RangeShiftR/man/figures/RS_repos.png)

RangeShifter is distributed with three user interfaces, each living in their own repo:

- the RangeShifter GUI (clickable Windows interface)*
- RangeShifter Batch Mode (command line interface)
- the RangeShiftR package (R interface)

All three share the same source code for the core simulation (i.e., the actual model), which lives in its own repo (RScore). Each of the interfaces keeps a copy of this core code in a subfolder called RScore, kept in sync with the RScore repo via a git subtree (see Git subtree usage section). 

*The RangeShifter GUI is currently being rewritten, and is not open source yet.


## Roles

#### Maintainers
- @JetteReeg
- @TheoPannetier

Maintainers are responsible for coordinating development efforts and ensuring that RangeShifter keeps building continuously.

#### Developers
Regular contributors and members of the RangeShifter development team

#### Contributors

Anyone who whishes to make changes to RangeShifter's code, including regular developers.

## Branching policy

This policy applies to RScore and all three RangeShifter interfaces.
RangeShifter uses the following branching structure:

- `main` is the default branch, where the stable releases live. Because it contains the version of RangeShifter that users normally interact with, it must be stable and build at all times.
Only maintainers should make changes to `main`, either directly for small changes (e.g. typo fixes), or by merging `develop` into `main` for any larger change. 
- `develop` is the development branch containing new features not yet made available to users.
Contributors are welcome to make changes to `develop`, but because this is the version that every contributor uses as a reference, one should ensure that new changes do not break `develop`.
If one happens to break `develop`, it should be their first priority to fix it.
For this reason, we recommend working from feature branches instead.
- Feature branches are created from `develop` by contributors to work on a new feature or other change, e.g. `cmake`, `mutualism`, etc. 
Contributors can also create their own branch, e.g. `theo` or `jette` to experiment with the code or implement miscellaneous changes.
Once a contributor deems their changes ready to be added to the development version, they should merge their changes from the feature branch into `develop`.
Optionally, we encourage contributors to seek a review from one or more developers and or maintainers by opening a pull request to merge their branch into develop.

### Contributing to RangeShifter core code

As mentioned above, the RangeShifter core code is located and maintained in thsi repo, which is currently only maintained by the RangeShifter development team. 
Any changes regarding the RangeShifter core code should be done in this repository and afterwards synced with all interfaces using the git subtree feature (see [Git subtree usage section](https://github.com/RangeShifter/RangeShiftR-package-dev/blob/development-guidelines/CONTRIBUTING.md#git-subtree-usage)). 

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
