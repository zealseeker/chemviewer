import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="chemviewer",
    version=__import__("chemviewer").__version__,
    author="Hongbin Yang",
    author_email="yanyanghong@163.com",
    description="A web based application to view 2D structures in a table",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zealseeker/chemviewer",
    packages=setuptools.find_packages(),
    entry_points={'console_scripts': ['chemviewer=chemviewer.run:main']},
    include_package_data=True,
    package_data = {
        'chemviewer': ['templates/*','static/*']
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    zip_safe=False,
)
