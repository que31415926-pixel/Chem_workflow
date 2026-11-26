# Software provenance

LICENSE applies to the original ChemFlow 2 engine, its new tests, examples and documentation in this directory. The MIT grant excludes vendor/, third-party dependencies, ORCA and research packs.

vendor/chemflow-v1/chemflow contains first-generation Python operators from this repository, based on f4a42aac9495ed40c411f73af98a53bbe55cf2b5, with compatibility and numerical fixes. That upstream revision has no standalone license file. Distribution of those operators outside this project requires the applicable permission of the original rights holders.

The private worker image installs scientific dependencies from conda-forge. The package metadata for the resolved versions carries the respective license texts. ORCA is supplied separately by the account holder under their agreement with the supplier. A private build uses the account holder's installation file. Neither its executable nor a prebuilt licensed image is distributed here.
