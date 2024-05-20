# Contributing

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways.

## Types of Contributions

### Report Bugs

Report bugs at <https://github.com/mn5hk/pyshbundle/issues>.

If you are reporting a bug, please include:

-   Your operating system name and version.
-   Any details about your local setup that might be helpful in troubleshooting.
-   Detailed steps to reproduce the bug.

### Found Bugs!!!

Look through the GitHub issues for bugs. Anything tagged with `bug` and
`help wanted` is open to whoever wants to implement it.

### Want New Functionality?

Look through the GitHub issues for features. Anything tagged with
`enhancement` and `help wanted` is open to whoever wants to implement it.

### Let's Improve Documentation

pyshbundle could always use more documentation,
whether as part of the official pyshbundle docs,
in docstrings, or even on the web in blog posts, articles, and such.<br>

Another prefered way is to create explainatory tutorials. Using the `PySHBundle` functions to explain the `Spherical Harmonics` and `GRACE Data Processing`.

### Feedback

The best way to send feedback is to file an issue at
<https://github.com/mn5hk/pyshbundle/issues>.

If you are proposing a feature:

-   Explain in detail how it would work.
-   Keep the scope as narrow as possible, to make it easier to implement.
-   Remember that this is a volunteer-driven project, and that contributions are welcome :)

## Get Started!

Ready to contribute? Here's how to set up pyshbundle for local development.

1.  Fork the pyshbundle repo on GitHub.

2.  Setup a seperate development environment.
    ```shell
    # clone the repo and fetch the dev branch
    $ git clone git@github.com:mn5hk/pyshbundle.git

    # creating a new virtual environment
    $ python3 -m venv <name-env>

    # install the dependencies from the requirements-dev file
    $ pip install -r ../pyshbundle/requirements-dev.txt

    # activate the virtual environment environment
    $ source </location-of-virt-env/name-env/bin/activate>
    ```

3. Build the latest repo in the development virtual environment.
    ```shell
    # install package into virtual environment
    $ pip install ../pyshbundle/dist/<required-version>.tar.gz

    # you also have the option to build the module using, be careful of 
    $ python setup.py sdist
    ```


4.  Create a branch for local development:

    ```shell
    $ git checkout -b name-of-your-bugfix-or-feature
    ```

    Now you can make your changes lo    cally.

5.  Commit your changes and push your branch to GitHub:

    ```shell
    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature
    ```

6.  Submit a pull request through the GitHub website.

## Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

TBD...
