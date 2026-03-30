# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Data management -- delegates to hitlist.

All data downloads, registration, and path resolution are handled by
hitlist. This module re-exports hitlist's API for backward compatibility.
"""

from __future__ import annotations

from hitlist.downloads import (
    FETCHABLE_DATASETS,
    MANUAL_DATASETS,
    available_datasets,
    data_dir,
    fetch,
    get_path,
    info,
    list_datasets,
    refresh,
    register,
    remove,
)

__all__ = [
    "FETCHABLE_DATASETS",
    "MANUAL_DATASETS",
    "available_datasets",
    "data_dir",
    "fetch",
    "get_path",
    "info",
    "list_datasets",
    "refresh",
    "register",
    "remove",
]
