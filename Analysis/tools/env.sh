#!/bin/bash

# Sample Tag and Probe environment variables

# get the project dir
tools_dir=$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd -P $tools_dir/.. && pwd)
export TNP=$project_dir
