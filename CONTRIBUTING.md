Contributing
============

Code
----

Please read the [standards document](doc/design/standards) before
contributing.

### Development

Make new branches for features `git checkout -b my_feature` and commit often
and push a little less often. Try to merge back to main branch as soon as
you have something that works.

#### Linting

For *C++* code make sure to use a static code checker like `cpplint` to
check the code for any style issues before submitting.

To install `cpplint`

```sh
# depending on your system one of these lines applies
pip install --user cpplint
pip install cpplint
python3 -m pip cpplint
python3 -m pip --user cpplint
```

### Commit Styling

The first line of the commit must be *at most* ~50 characters long and
should start with either.

- `Feat:` For new feature.
- `Fix:` For bug fix.
- `Merge:` For merging.

The commit first line must be in *present* tense so that anyone picking a
commit hash can easily read what they are enabling.

### Pull Requests

Make sure you have linted and checked your code before asking for a pull
request. Another member must check the code and approve it before merge.

Issues
------

*Issues* are reporting bugs, feature requests or goals for the project. In
order to submit an issue make sure it follows the [issue
template](.github/ISSUE_TEMPLATE).
