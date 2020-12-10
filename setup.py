from setuptools import setup, find_packages

setup(
    name="recover_data_lake",
    description='The ingestion pipelines for the RECOVER data lake',
    packages=find_packages(include=["recover_data_lake", "recover_data_lake.*"]),
    version="0.1")
