.. _developer_guide:

Developer guide
=======================

This section is under construction.


..
   Thank you for considering contributing to this project! Here are some guidelines to help you get started:

   Developer Installation Guide
   ----------------------------

   This installation guide is intended for developers who wish to contribute to or modify the barotropy source code. It assumes that the developer is using a Linux distribution or Windows with Git Bash terminal to have access to Git and Linux-like commands.

   1. **Fork the repository:**

      - Navigate to the `project's GitHub page <https://github.com/turbo-sim/barotropy>`_.
      - Click the "Fork" button in the upper right corner of the repository page to create a copy of the repository under your own GitHub account.


   2. **Clone the forked repository:**

      - Open your terminal.
      - Run the following command, replacing `<your-username>` with your GitHub username:

      .. code-block:: bash

         git clone https://github.com/<your-username>/<repository-name>.git

      - Navigate into the cloned repository:

      .. code-block:: bash

         cd <repository-name>

   3. **Create a dedicated Conda virtual environment for barotropy development**:

      - Check that conda is installed:

      .. code-block:: bash

         conda list

      - If not conda is installed, `install conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_.
      - Create dedicated virtual environment for barotropy package:

      .. code-block:: bash

         conda env create --name barotropy_env python=3.11

   4. **Activate the newly created virtual environment**:

      .. code-block:: bash

         conda activate barotropy_env

   5. **Install Poetry to manage dependencies**:

      .. code-block:: bash

         conda install poetry

      Poetry is a powerful dependency manager that offers separation of user and developer dependencies, ensuring that only the necessary packages are installed based on the user's intent. Additionally, it simplifies the process of adding, updating, and removing dependencies, making it easier to maintain the project's requirements.

   6. **Use Poetry to install the required dependencies for barotropy development**:

      .. code-block:: bash

         poetry install

   7. **Verify the installation by running the following command**:

      .. code-block:: bash

         python -c "import barotropy; barotropy.print_package_info()"

      If the installation was successful, you should see the barotropy banner and package information displayed in the console output.


   Pull request guidelines
   -------------------------

   Please follow these steps to submit a pull request.

   1. **Create a branch in your forked repository**:

      - Open your terminal in the projects root.
      - Create branch:

      .. code-block:: bash

         git checkout -b <feature-name>

   2. **Make your changes**:

      - Implement your feature or bugfix.


   3. **Commit your changes**:

      .. code-block:: bash 

         git commit -m "Description of changes"

   4. **Push to your fork**: 

      .. code-block:: bash

         git push origin feature-name

   5. **Open a pull request**: 

      - Go to your fork on GitHub and click the "New pull request" button.



   Reporting issue
   ----------------

   If you find a bug or have a feature request, please open an issue in the Github project page and follow the provided templates.

   CI/CD Pipeline
   --------------

   barotropy uses GitHub Actions to automate its Continuous Integration and Continuous Deployment (CI/CD) processes.

   Automated Testing
   ^^^^^^^^^^^^^^^^^

   The ``ci.yml`` action is triggered whenever a commit is pushed to the repository. This action runs the test suite on both Windows and Linux environments, ensuring the code's compatibility and correctness across different platforms.

   Package Publishing
   ^^^^^^^^^^^^^^^^^^

   barotropy utilizes the ``bumpversion`` package to manage versioning and release control. To increment the version number, use the following command:

   .. code-block:: bash

      bumpversion patch  # or minor, major

   After bumping the version, push the changes to the remote repository along with tags to signify the new version:

   .. code-block:: bash

      git push origin --tags

   If the tests pass successfully, the package is automatically published to the Python Package Index (PyPI), making it readily available for users to install and use.

   Documentation Deployment
   ^^^^^^^^^^^^^^^^^^^^^^^^

   barotropy automates the deployment of documentation using the ``deploy_docs`` action. This action builds the Sphinx documentation of the project and publishes the HTML files to GitHub Pages each time that a new commit is pushed to the remote repository. By automating this process, barotropy ensures that the project's documentation remains up-to-date and easily accessible to users and contributors.