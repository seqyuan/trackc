"""Utility functions and classes
"""
import sys
import inspect
import warnings
import pandas as pd

class GenomeRegion:
    region = None
    chrom = None
    start = None
    end = None
    length = None
    isReverse = False
    
    def __init__(self, region: str):
        self.region = region
        tmp = region.split(":")
        self.chrom = tmp[0]
        if len(tmp)==2:
            self.start = int(tmp[1].split("-")[0])
            self.end = int(tmp[1].split("-")[1])
            if self.start > self.end:
                self.isReverse = True
            self.length = abs(self.start - self.end)
        
    def GenomeRegion2df(self):
        region4coolFetch = self.chrom + ":" + str(self.start) + '-' + str(self.end)
        if self.start == None:
            region4coolFetch = self.chrom
        else:
            if self.start > self.end:
                region4coolFetch = self.chrom + ":" + str(self.end) + '-' + str(self.start)

        df = pd.DataFrame(
            {'chrom':[self.chrom], 
            'start':[self.start],
            'end':[self.end],
            'isReverse':[self.isReverse],
            'region4coolFetch': [region4coolFetch]
            }, 
            index=[self.region])

        return df

my23colors = ['#53868B','#00F5FF','#C1FFC1','#0000FF','#7B68EE',
                  '#CDCD00','#FFF68F','#CD9B1D','#8B658B','#FF6A6A','#8B3A3A',
                  '#1E90FF','#FF69B4','#8DB6CD','#CAE1FF','#EECFA1','#8B7B8B',
                  '#4F4F4F','#FF4500','#BC8F8F','#FFA500','#228B22','#8B4513']


def getdoc(c_or_f: Union[Callable, type]) -> Optional[str]:
    if getattr(c_or_f, '__doc__', None) is None:
        return None
    doc = inspect.getdoc(c_or_f)
    if isinstance(c_or_f, type) and hasattr(c_or_f, '__init__'):
        sig = inspect.signature(c_or_f.__init__)
    else:
        sig = inspect.signature(c_or_f)

    def type_doc(name: str):
        param: inspect.Parameter = sig.parameters[name]
        cls = getattr(param.annotation, '__qualname__', repr(param.annotation))
        if param.default is not param.empty:
            return f'{cls}, optional (default: {param.default!r})'
        else:
            return cls

    return '\n'.join(
        f'{line} : {type_doc(line)}' if line.strip() in sig.parameters else line
        for line in doc.split('\n')
    )


def _one_of_ours(obj, root: str):
    return (
        hasattr(obj, "__name__")
        and not obj.__name__.split(".")[-1].startswith("_")
        and getattr(
            obj, '__module__', getattr(obj, '__qualname__', obj.__name__)
        ).startswith(root)
    )

def descend_classes_and_funcs(mod: ModuleType, root: str, encountered=None):
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
            yield from descend_classes_and_funcs(obj, root, encountered)


def annotate_doc_types(mod: ModuleType, root: str):
    for c_or_f in descend_classes_and_funcs(mod, root):
        c_or_f.getdoc = partial(getdoc, c_or_f)


def _doc_params(**kwds):
    """\
    Docstrings should start with "\" in the first line for proper formatting.
    """

    def dec(obj):
        obj.__orig_doc__ = obj.__doc__
        obj.__doc__ = dedent(obj.__doc__).format_map(kwds)
        return obj

    return dec


def _check_array_function_arguments(**kwargs):
    """Checks for invalid arguments when an array is passed.

    Helper for functions that work on either AnnData objects or array-likes.
    """
    # TODO: Figure out a better solution for documenting dispatched functions
    invalid_args = [k for k, v in kwargs.items() if v is not None]
    if len(invalid_args) > 0:
        raise TypeError(
            f"Arguments {invalid_args} are only valid if an AnnData object is passed."
        )


# --------------------------------------------------------------------------------
# xx stuff
# --------------------------------------------------------------------------------


