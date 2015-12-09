#!/usr/bin/env python
"""Classes to assist with GFF parsing. It tries to keep track of the partent/
child relationships and allows for fetching of parent or children features.

Contains the following classes:
    GFFError:   An error when a queried ID is not found in the GFF data.
    GFFFeature: Reads a line from a GFF file and parses out the feature info
    GFFHandler: Handles the parsing of the whole GFF file and allows access
                to parents, children, or "siblings" of individual features.
"""


class GFFError(LookupError):
    """Raised when a user asks for a feature ID that is not present in the GFF
    data."""


class GFFFeature(object):
    """This reads the data out of a single line of a GFF. It is a *very* ugly
    class, but GFF is a *very* ugly file format, and I don't know of a better
    way to handle this data."""
    __slots__ = [
        'seqid',
        'source',
        'type',
        'start',
        'end',
        'score',
        'strand',
        'phase',
        'ID',
        'Name',
        'Alias',
        'Parent',
        'Target',
        'Gap',
        'Derives_from',
        'Note',
        'Dbxref',
        'Ontology_term',
        'Is_circular'
        ]

    def __init__(self, line):
        self.seqid = None
        self.source = None
        self.type = None
        self.start = None
        self.end = None
        self.score = None
        self.strand = None
        self.phase = None
        self.ID = None
        self.Name = None
        self.Alias = None
        self.Parent = None
        self.Target = None
        self.Gap = None
        self.Derives_from = None
        self.Note = None
        self.Dbxref = None
        self.Ontology_term = None
        self.Is_circular = None
        self.read_gff_line(line)

    def read_gff_line(self, line):
        fields = line.strip().split('\t')
        #   Fixed fields are the first eight
        self.seqid = fields[0]
        self.source = fields[1]
        self.type = fields[2]
        self.start = int(fields[3])
        self.end = int(fields[4])
        self.score = fields[5]
        self.strand = fields[6]
        self.phase = fields[7]
        #   For the next field, we split it and look for keywords
        attrib = fields[8].split(';')
        for a in attrib:
            if a.startswith('ID'):
                self.ID = a.split('=')[1]
            elif a.startswith('Name'):
                self.Name = a.split('=')[1]
            elif a.startswith('Alias'):
                alias = a.split('=')[1]
                self.Alias = alias.split(',')
            elif a.startswith('Parent'):
                parent = a.split('=')[1]
                self.Parent = parent.split(',')
            elif a.startswith('Target'):
                self.Target = a.split('=')[1]
            elif a.startswith('Gap'):
                self.Gap = a.split('=')[1]
            elif a.startswith('Derives_from'):
                self.Derives_from = a.split('=')[1]
            elif a.startswith('Note'):
                note = a.split('=')[1]
                self.Note = note.split(',')
            elif a.startswith('Dbxref'):
                dbxref = a.split('=')[1]
                self.Dbxref = dbxref.split(',')
            elif a.startswith('Ontology_term'):
                ont = a.split('=')[1]
                self.Ontology_term = ont.split(',')
            elif a.startswith('Is_circular'):
                self.Is_circular = a.split('=')[1]


class GFFHandler(object):
    """Handles the parsing of the entire GFF file, line by line. Defines
    methods that allow for retrieval of specific features, or groups of
    features.

    Contains the following methods:
        get_feature(self, feat_id)
            Retuns a single feature, or raises an exception when the ID is not
            found

        get_parents(self, feat_id, type=None)
            Returns the parents of a feature, or raises an exception when the
            ID is not found. If the type keyword is specified, then only
            features of a certain type are returned. If there are no features
            of the supplied type, then empty list is returned. Defaults to
            returning all parents. This will only traverse up one level.

        get_children(self, feat_id, type=None)
            Returns the children of a feature, or raises an exception when the
            ID is not found. If the type keyword is specified, then only child
            features of the specified type are returned. If there are no child
            features of that type, then empty list is returned. Defaults to
            returning all children. This will only traverse down one level.

        get_siblings(self, feat_id, type=None)
            Returns the "siblings" (features that share parents), or raises an
            exception when the ID is not found. If the type keyword is
            specified, then only features of that type are returned. If there
            are no sibling features of that type, then empty list is returned.
            Defaults to returning all siblings. This does not check parents
            of parents or children of children.

        chrom_features(self, chrom, left, right, type=None)
            Returns all features on the specified chromosome, or raises an
            exception when that chromosome is not found. If the type keyword is
            specified, then only features of that type are returned. If the
            left keyowrd is specified, then features that have bases to the
            _right_ of the specified base are returned. If the right keyword is
            specified, then features with bases to the _left_ of the offset
            are returned. These effctively limit the search to a certain
            region. If there are no features of the specified type in the given
            region, then empty list is returned. Defaults to all types and the
            entire chromosome.
    """
    gff_data = []
    gff_ids = []

    def __init__(self):
        pass

    def gff_parse(self, fname):
        with open(fname, 'r') as f:
            for line in f:
                #   Skip lines starting with #; they are directives or comments
                if line.startswith('#'):
                    continue
                #   Skip blank lines. I guess we should support CRLF lines too
                elif line in ('\n', '\r\n'):
                    continue
                else:
                    g = GFFFeature(line)
                    self.gff_data.append(g)
                    self.gff_ids.append(g.ID)

    def get_feature(self, feat_id):
        if feat_id not in self.gff_ids:
            raise GFFError('ID {fid} not found in GFF!'.format(fid=feat_id))
        else:
            for f in self.gff_data:
                if f.ID == feat_id:
                    return(f)

    def get_parents(self, feat_id, feat_type=None):
        if feat_id not in self.gff_ids:
            raise GFFError('ID {fid} not found in GFF!'.format(fid=feat_id))
        else:
            parents = []
            targeted_features = []
            #   We search through the list and save features that match the
            #   query ID
            for f in self.gff_data:
                if f.Parent and feat_id == f.ID:
                    if feat_type:
                        if f.type == feat_type:
                            parents += f.Parent
                    else:
                        parents += f.Parent
            #   We then search AGAIN to get the features that match the parent
            #   feature ID. This is not very efficient, but given the structure
            #   of the data, it works.
            for f in self.gff_data:
                if f.ID in parents:
                    targeted_features.append(f)
            return(targeted_features)

    def get_children(self, feat_id, feat_type=None):
        if feat_id not in self.gff_ids:
            raise GFFError('ID {fid} not found in GFF!'.format(fid=feat_id))
        else:
            targeted_features = []
            #   Searching for child features is easy, since we just have to see
            #   if our query ID is in the Parent GFF attributes.
            for f in self.gff_data:
                if f.Parent and feat_id in f.Parent:
                    if feat_type:
                        if f.type == feat_type:
                            targeted_features.append(f)
                    else:
                        targeted_features.append(f)
            return(targeted_features)

    def get_siblings(self, feat_id, feat_type=None):
        if feat_id not in self.gff_ids:
            raise GFFError('ID {fid} not found in GFF!'.format(fid=feat_id))
        #   This is pretty easy, too. We just take the query ID, find its
        #   parents, and then find all children of those parents. These are the
        #   "siblings"
        parents = self.get_parents(feat_id, feat_type=feat_type)
        sibs = []
        for p in parents:
            sibs += self.get_children(p.ID, feat_type=feat_type)
        return(sibs)

    def chrom_features(self, chrom, left=None, right=None, feat_type=None):
        targeted_features = []
        filtered_features = []
        for f in self.gff_data:
            if f.seqid == chrom:
                if feat_type:
                    if f.type == feat_type:
                        targeted_features.append(f)
                else:
                    targeted_features.append(f)
        #   Now we go through the list of targeted features and filter based
        #   on position
        if left or right:
            #   if left is not specified, set it to 0
            if not left:
                left = 0
            #   Same with right, except make it huge
            if not right:
                right = 1e99
            left = int(left)
            right = int(right)
            for t in targeted_features:
                if t.start in range(left, right) or t.end in range(left, right):
                    filtered_features.append(t)
        else:
            filtered_features = targeted_features
        return(filtered_features)

    def overlapping_feature(self, chrom, pos, feat_type=None):
        targeted_features = []
        for f in self.gff_data:
            if f.seqid == chrom:
                if int(pos) in range(f.start, f.end):
                    if feat_type:
                        if f.type == feat_type:
                            targeted_features.append(f)
                    else:
                        targeted_features.append(f)
        return targeted_features
