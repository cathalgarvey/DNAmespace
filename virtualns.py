'''virtualns - A Virtual Namespace object providing a dict-like interface.
by Cathal Garvey
Licensed under the GNU Affero GPL, the license text of which can be accessed
as virtualns.license in a python prompt or script.
Designed as part of the DNAmespace project.
'''
import keyword
from gnulicenses import Affero as license

class nsdict(object):
    '''Mimics a dict but has no exposed methods, and exposes all keys as attributes.
    This allows use of this object as a "Virtual Namespace" with a dict-like interface.
    Keys that are added to this object using the "dict[key] = value" idiom can be
    accessed as dict.key, and iPython autocompletion will work on these keys.

    Because there are more pitfalls to using attributes versus dict keys, such
    as bugs due to setting a reserved keyword as an attribute, this object takes
    measures to try and prevent this by suffixing disallowed or potentially
    dangerous attribute names with underscores. Some of this magic happens
    in the __init__ method, so subclasses of this object should call the parent
    __init__ if they define their own __init__ method.

    Barring a few kludges to make things work as expected, this object
    has all of the expected dict methods, but prepended with an underscore to "hide"
    them from attribute listing and autocompletion in iPython etc.

    For example, this object provides a "_keys()" attribute, which returns
    a list of keys (an actual list and not a dict_keys, due to kludges),
    and a pop method, which works as normal.

    In most cases, these underscored dict methods are simply passed with
    all arguments to the object's underlying __dict__ attribute.

    The only significant difference in behaviour between this class and a dict
    is that it cannot be directly passed to another dict's "update" method, as
    its equivalent of self.keys() is obfuscated. Instead, a compatibility method
    is added: self._asupdate(), which returns a normal dict containing all
    significant key:value pairs from self.__dict__ (that is, those that are not
    in self._method_keywords).'''

    def __init__(self, somedict=None, autofix=False):
        '''If autofix is set to true, then conflicting keywords are
        automatically suffixed with underscores until they no longer conflict.'''
        self._autofix = bool(autofix)
        # Need a list of methods so we can rename keys that conflict.
        # Keep this call second last in __init__, right before dict import.
        self._method_keywords = dir(self)
        self._method_keywords.append("_method_keywords")
        # If provided with a dict, then use self._update to import it.
        if somedict:
            # Could be strict on type input, but that would prevent
            # passing of dict subclasses and iterables.
            # Let __dict__ raise whatever TypeError may result in case
            # of non-updateable input.
            self._update(somedict)

    def __test_conflict__(self, key):
        'Tests a new key for conflict with reserved keywords or required object methods/attributes.'
        return (key in keyword.kwlist) or (key in self._method_keywords)

    def __setitem__(self, key, value):
        'Set properties for all new dict keys. Reserved words are suffixed with "_".'
        # Make sure the property-ified key won't collide with reserved:
        if self.__test_conflict__(key):
            if self._autofix:
                while self.__test_conflict__(key):
                    # If the key conflicts, keep appending underscores.
                        key = key+"_"
            else:
                raise KeyError(("Key to be set conflicts with a reserved "
                    "keywords or required object attribute."))
        self.__dict__[key] = value

    def __getitem__(self, key):
        'Should behave like a similar call to dict.'
        return self.__dict__[key]

    def __delitem__(self, key):
        'Should behave like a similar call to dict.'
        del(self.__dict__[key])

    def _clear(self):
        return self.__dict__.clear()

    def _copy(self):
        return self.__dict__.copy()

    def _fromkeys(self, *args, **nargs):
        return self.__dict__(*args, **nargs)

    def _get(self, *args, **nargs):
        return self.__dict__.get(*args, **nargs)

    def _items(self):
        return self.__dict__.items()

    def _keys(self):
        'Technically this returns a list, not a dict_keys object.'
        keys = list(self.__dict__.keys())
        for item in self._method_keywords:
            try:
                del(keys[keys.index(item)])
            except ValueError: # When something isn't in the list.
                pass
        return keys

    def _mro(self):
        return self.__dict__.mro()

    def _pop(self, *args, **nargs):
        return self.__dict__.pop(*args, **nargs)

    def _popitem(self, *args, **nargs):
        return self.__dict__.popitem(*args, **nargs)

    def _setdefault(self, *args, **nargs):
        return self.__dict__.pop(*args, **nargs)

    def _update(self, *args, **nargs):
        # nsdicts can be updated with other dicts as normal. However, to
        # perform the reverse operation, use dict.update(nsdict._asupdate()).
        return self.__dict__.update(*args, **nargs)

    def _values(self):
        return self.__dict__.values()

    def _asupdate(self):
        '''Allows passing of an nsdict object to the update method of a normal
        dict, by returning a redacted copy of the nsdict's __dict__ attribute,
        removing any key/value pairs that match the self._method_keywords list.'''
        compatible_dict = self._copy()
        # Clear out methods and attribute names, spare only namespace items.
        for key in self._method_keywords:
            try:
                del(compatible_dict[key])
            except KeyError:
                pass
        return compatible_dict
