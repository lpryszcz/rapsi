
class htmlTable: 
    """General Table formatter. It can be used to create a grid structure
    using python syntaxis, and convert it into HTML format.
    by jhcepas@gmail.com
    ##
    # Example: 

    T = htmlTable()
    T.add_cell(0,"Phylome")
    T.add_cell(0,"Description")
    T.add_cell(1,"Human Phylome")
    T.add_cell(1,"2007 Genome Biology paper.")

    # Set CSS style used for different rows in table

    T.table_style = "gtable"
    T.header = True 

    print T.asHTML()
    """
    def __init__(self):
        self.rows = []
        self.header = True
        self.table_style = ""
        self.header_style = ""
        self.content_style = ""
        self.td_flags = ""
        self.row_classes = {}
        self.cell_classes= {}

    def add_cell(self, row_index, content, classname=None):
        """Add the given object (any text or HTML code) as the content of
        a cell in the row specified by the "row_index" argument. Cells
        are arranged in a stack-like way. """
        # Creates the grid until fitting the given index
        while len(self.rows)<=row_index:
            self.rows.append([])
        self.rows[row_index].append(str(content))
        #if classname is not None: 
        #    self.cell_classes[row_index][len(self.rows[row_index]-1)] = classname
            
    def remove_column(self, column_index):
        """Remove column from table"""
        for row_index, row in enumerate(self.rows):
            if len(row)<=column_index:
                continue
            #remove column content
            row.pop(column_index)
            self.rows[row_index] = row
            #drop formatting classes
            #if self.cell_classes[row_index]:

    def asHTML(self):
        """Returns the HTML representation of the table object."""
        html = ['<table align=center class="%s">' % self.table_style]
        for r in self.rows:
            html.append(' <tr>\n  ')
            for cell in r:
                if self.rows.index(r) == 0 and self.header:
                    html[-1] += '<th %s >%s</th>' %(self.td_flags, cell)
                else:
                    html[-1] += '<td %s >%s</td>' %(self.td_flags, cell)
            html.append(" </tr>")
        html.append("</table>")
        return "\n".join(html) + "\n"

    def asTXT(self):
        """Returns the TEXT representation of the table object."""
        txt = []
        for r in self.rows:
            if self.rows.index(r) == 0 and self.header:
                txt.append('#'+'\t'.join(r))
            else:
                txt.append('\t'.join(r))
        return "\n".join(txt) + "\n"
