# wick

`wick` is an R package to access and query the `CANDEL` datastore.

`wick` communicates with `CANDEL` using [datalogr](https://github.com/CANDLEbio/datalogr). Please refer to this README for information on database access and the DSL that enables you to write Datalog queries from R

`wick` requires a CANDEL instance. An example setup can be found at the [datalog-json-parser](https://github.com/CANDELbio/datalog-json-parser/tree/master/examples/client-server)


## Pre-defined queries

`wick` contains a number of pre-defined queries. For each query, the function `get_QUERYNAME` runs the query against the database and returns the query result, while the function `QUERYNAME` returns the data representing the query itself


Copyright Parker Institute for Cancer Immunotherapy, 2022

Licensing :
This script is released under the Apache 2.0 License
Note however that the programs it calls may be subject to different licenses.
Users are responsible for checking that they are authorized to run all programs
before running this script.