from skbuild import setup

setup(
    packages=["geomat"],
    package_dir={"": "src"},
    cmake_install_dir="src/geomat",
    include_package_data=True,
)
