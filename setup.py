from setuptools import setup, find_packages

setup(
    name="genofunc",
    version="0.1.0",
    packages=find_packages(),
    url="https://github.com/xiaoyu518/genofunc",
    license="MIT",
    entry_points={"console_scripts": ["genofunc = genofunc.__main__:main"]},
    test_suite="nose.collector",
    tests_require=["nose >= 1.3"],
    install_requires=[
        "biopython>=1.70",
        "numpy>=1.18",
        "pandas>=0.24.2",
        "nextstrain-augur>=13.1.2",
        "mappy>=2.24"
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
