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

⚠️ If you wish to propose a change to the core code of the simulation, please do so *in the [RScore](https://github.com/RangeShifter/RScore) repo*, rather than in the RScore folder of either interface.

*The RangeShifter GUI is currently being rewritten, and is not open source yet.

## RangeShiftR package

### Issues

Issues should be used for reporting technical problems with the R package or suggest improvements e.g. in the documentation. 
Questions regarding the application and settings of the RangeShiftR should be asked in the [discussion forum](https://github.com/RangeShifter/RangeshiftR-tutorials/discussions).

#### Create a new issue

If you encounter a technical problem in the RangeShiftR package or would like to suggest an improvement to the package, please search if [a related issue already exists](https://github.com/RangeShifter/RangeShiftR-package-dev/issues). If you can't find a related issue, you can open a new issue using a relevant [issue form](https://github.com/RangeShifter/RangeShiftR-package-dev/issues/new/choose).

#### Solve an issue

As a general rule, the RangeShiftR development team will take care of opened issues. However, feel free to scan through [existing issues](https://github.com/RangeShifter/RangeShiftR-package-dev/issues) and leave comments. 

### Fixing a bug

If you found a bug in the RangeShiftR package files, which isn't already reported in an [existing issue](https://github.com/RangeShifter/RangeShiftR-package-dev/issues), please report the bug using the dedicated [issue form](https://github.com/RangeShifter/RangeShiftR-package-dev/issues/new/choose).
If you were able to solve the bug already, you can suggest the fix in the issue discussion and create a pull request. 
The pull request will then be reviewed by the RangeShiftR development team and eventually be merged.

### Adding a new feature

Any suggestions how to improve the R package? Or do you think, we missed a feature in the R package interface? Please check if someone else [already suggested the feature](https://github.com/RangeShifter/RangeShiftR-package-dev/issues) and use the [dedicated form](https://github.com/RangeShifter/RangeShiftR-package-dev/issues/new/choose) to suggest a new feature. 
You can also contact the RangeShiftR development team directly (rangeshiftr@uni-potsdam.de).
