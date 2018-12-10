from setuptools import setup

if __name__ == "__main__":

    setup(
        name="kaplan",
        version="0.0.0",
        description="Conformer identification package.",
        url="https://github.com/PeaWagon/Kaplan",
        author="Jen Garner",
        author_email="garnej2@mcmaster.ca",
        package_dir={"kaplan": "kaplan"},
        requires=["numpy"],
        )
