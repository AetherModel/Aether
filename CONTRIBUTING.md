# Contributing

Please read the [standards document](doc/design/standards) before
contributing.

- [Contributing](#contributing)
  - [Code](#code)
    - [Development](#development)
      - [Using AStyle](#using-astyle)
      - [Linting](#linting)
  - [Commit Styling](#commit-styling)
  - [Pull Requests](#pull-requests)
  - [Issues](#issues)

---

## Code

### Development

Make new branches for features `git checkout -b my_feature` and commit often
and push a little less often. Try to merge back to main branch as soon as
you have something that works.

#### Using AStyle

We have started using astyle to format the code.  Please see
<http://astyle.sourceforge.net/>. There is a style file in the root directory of
Aether, called .astylerc.  To run this, do:

AStyle --project=.astylerc src/*.cpp

on WSL with Ubuntu, the command seems to be:

astyle --options=.astylerc src/*.cpp

Astyle can also be used within VSCode as an extension. Install the "Astyle" extension (author: Chieh Yu).  Then go to settings (the "Manage" wheel on the bottom left, then click on "settings"). Under "Extensions", click on "Astyle Configuration".  Within that, you can set the location of the astylerc file. Click on "Edit in settings.json". Modify the file so it has the following line:
    "[cpp]": {
        "editor.defaultFormatter": "chiehyu.vscode-astyle",
        "editor.formatOnSave": true
    },
    "astyle.astylerc": "${workspaceRoot}/.astylerc"

This should then allow you to format the code using the defaults that we have set in the astyle file. This can be done with (Ctrl-Shift-I for Linux, Alt-Shift-F for Windows) and upon saving the file - this reformats the whole file to put it into the astyle format.

#### Linting

For *Python*,
[flake8](https://flake8.pycqa.org/en/latest/) is a good option.  This may by installed using pip.


## Commit Styling

The first line of the commit must be *at most* ~50 characters long and
should start with either.

- `FEAT:` For new feature.
- `BUG:` For bug fix.
- `MERGE:` For merging.
- `DOC:` For documentation update.
- `TEST:` For the addition or modification of tests.
- `STY:` For a style update (e.g., linting).
- `DEPREC:` Deprecate something, or removee a deprecated object.
- `REVERT:` Revert an earlier commit.
- `MAINT:` For maintenance such as refactoring, typos, etc.

The commit first line must be in *present* tense so that anyone
picking a commit hash can easily read what they are enabling. For more
information check out [conventional commit
messages](https://www.conventionalcommits.org/en/v1.0.0/).

For example,

*do:*

```github
FEAT: Hydrostatic density implementation.
```

*don't:*

```github
Implemented hydrostatic density. (feature)
```

## Pull Requests

Make sure you have linted and checked your code before asking for a
pull request. Before requesting a review, ensure the pull request
check list has been completed.  Another member must check the code and
approve it before merge.

## Issues

*Issues* are reporting bugs, feature requests, or goals for the
project. In order to submit an issue make sure it follows the [issue
template](.github/ISSUE_TEMPLATE).  Please search through the existing
issues before submitting a new one, as someone else may have already
come accross and reported the problem you've encountered.
