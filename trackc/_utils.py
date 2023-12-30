"""Utility functions and classes"""
from functools import partial, wraps
from types import MethodType, ModuleType
from typing import Any, Callable, Dict, Mapping, Optional, Tuple, Union
from weakref import WeakSet

import pandas as pd


def _one_of_ours(obj, root: str):
    return (
        hasattr(obj, "__name__")
        and not obj.__name__.split(".")[-1].startswith("_")
        and getattr(
            obj, "__module__", getattr(obj, "__qualname__", obj.__name__)
        ).startswith(root)
    )


def _descend_classes_and_funcs(mod: ModuleType, root: str, encountered=None):
    if encountered is None:
        encountered = WeakSet()
    for obj in vars(mod).values():
        if not _one_of_ours(obj, root):
            continue
        if callable(obj) and not isinstance(obj, MethodType):
            yield obj
            if isinstance(obj, type):
                for m in vars(obj).values():
                    if callable(m) and _one_of_ours(m, root):
                        yield m
        elif isinstance(obj, ModuleType) and obj not in encountered:
            encountered.add(obj)
            yield from _descend_classes_and_funcs(obj, root, encountered)
