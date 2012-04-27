import re

# BRETT: Minor mods here made to update it
#   re NOT pre
#   __doc__ not required

"""
This is a simple HTML table parser, modelled loosely after Perl's HTML::TableExtract.

The parse() function accepts an HTML function and returns a list of lists.

See http://bebop.bigasterisk.com/python for more Python stuff.
This file is available under the GPLv2 license.
"""

__license__ = 'GPLv2'
__author__ = 'David McClosky (dmcc@bigasterisk.com)'
__version__ = '2.0'
__all__ = ['clean', 'parse']

# TODO: distutils-fication

def tag_maker(tag, close_tag=None):
    "Make regexes that match HTML tags"
    if not close_tag: close_tag = tag
    return re.compile(r"<%s[^>]*>(.*?)</%s[^>]*>" % (tag, close_tag), 
                      re.I | re.M | re.S)

def ent_maker(ent):
    "Make regexes that match entities"
    return re.compile("&%s;?" % ent, re.I)

tbl = tag_maker('table')
hdr = tag_maker('th', 't(?:r|h|able)')
row = tag_maker('tr', 't(?:r|able)')
cel = tag_maker('td', 't(?:d|r|able)')
tag = re.compile(r"<[^>]*>", re.M | re.S)

ent = {}

for (old, new) in (('amp', '&'), ('lt', '<'), ('gt', '>'), ('nbsp', ' ')):
    ent[ent_maker(old)] = new

def clean(str):
    'Removes tags and surrounding whitespace'
    t = tag.sub('', str)
    for pat, repl in ent.items():
        t = re.sub(pat, repl, t)
    return t.strip()

def parse(str, head=None, cleaner=clean):
    """Parse a string, given optionally a list of header entries:

    parse('<HTML>...(snip)...</HTML>', ['first name', 'last name'])
    
    A function to cleanup the HTML may be specified if you don't like
    the behavior of the original one.  It will return a list of lists."""
    newtbl = []
    tbls = tbl.findall(str)
    for t in tbls:
        header = [cleaner(h) for h in hdr.findall(t)]
        if head and head != header: continue
        else:
            newtbl.append(header)
            rows = row.findall(str)
            for r in rows:
                newrow = []
                cells = cel.findall(r)
                for c in cells: newrow.append(clean(c))
                newtbl.append(newrow)
            return newtbl
    else: return None


