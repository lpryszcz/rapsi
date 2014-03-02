
import sqlite3

class Taxonomy(object):
    """Handle NCBI taxonomy stored in sqlite3 db."""
    def __init__(self, db):
        self.cnx = sqlite3.connect(db)
        self.cur = self.cnx.cursor()
        
    def __getitem__(self, taxid):
        return self.get_taxa_info(taxid)

    def get_taxa_info(self, taxid):
        """Return taxa info: parent taxid, species name and rank."""
        cmd = 'select parent_id, name, rank from taxa_info where taxid=%s'%taxid
        self.cur.execute(cmd)
        result = self.cur.fetchone()
        if result:
            return result
        else:
            return None, None, None

    def gi2taxid(self, gi):
        """Fetch taxid from gi2taxid table."""
        self.cur.execute('select taxid from gi2taxid where gi=%s'%gi)
        result = self.cur.fetchone()
        if result:
            return result[0]
        else:
            return None
            
    def taxid2lineage(self, taxid, maxRank="no rank"):
        """Return taxid linneage. Stop after reaching maxRank (root)."""
        parents = []
        rank = ""
        while taxid!=1 and rank!=maxRank:
            taxid, name, rank = self.get_taxa_info(taxid)
            if not taxid:
                break
            parents.append(taxid)
        return parents

    def collapse_taxa(self, taxids, maxRank="family"):
        """Return collapsed taxa on the level of species or genus."""
        #get parents
        taxid2childs = {}
        taxid2irank  = {}
        #convert to set
        if not type(taxids) is set:
            taxids = set(taxids)
        for taxid in taxids:
            #get species, species group and genus for each taxa
            #and generate childs
            for i, p in enumerate(self.taxid2lineage(taxid, maxRank)):
                if p not in taxid2childs:
                    taxid2childs[p] = []
                    #store rank as int
                    taxid2irank[p]  = i
                #store child
                taxid2childs[p].append(taxid)
        #get common parent, if not return None (likely low-complexity)
        ptaxid = None
        irank  = 99
        for p in sorted(taxid2childs, key=lambda k: len(taxid2childs[k]), reverse=True):
            #skip parentals not covering all childs
            if taxids.difference(taxid2childs[p]):
                continue
            #replace ptaxid only if deeper rank (ie species instead of genus)
            if taxid2irank[p] < irank:
                ptaxid, irank = p, taxid2irank[p]
        return ptaxid
