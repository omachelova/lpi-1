import cnf

from typing import List, Mapping
Valuation = Mapping[str, bool]

#Literal
# Clause
#  Cnf
#   Formula
#    BinaryFormula
#     Implication
#     Equivalence

class Formula(object):
    def __init__(self, subs : List['Formula']= []) -> None:
        self.m_subf = subs # type: List[Formula]
    def subf(self) -> List['Formula']:
        return self.m_subf
    def isSatisfied(self, v : Mapping[str,bool]) -> bool:
        return False
    def toString(self) -> str:
        return "INVALID"
    def __str__(self) -> str:
        return self.toString()
    def __repr__(self) -> str:
        return self.__class__.__name__ + '(' + ','.join([ repr(f) for f in self.subf()]) + ')'

#    def toCnf(self) -> cnf.Cnf:
#        """ Vrati reprezentaciu formuly v CNF tvare. """
#        # TODO
#        return cnf.Cnf()

class Variable(Formula):
    def __init__(self, name : str) -> None:
        Formula.__init__(self)
        self._name = name
    def name(self) -> str:
        return self._name
    def isSatisfied(self, v : Mapping[str, bool]) -> bool:
        return v[self.name()]
    def toString(self) -> str:
        return self.name()
    def __repr__(self) -> str:
        return "Variable(%r)" % (self.name(),)
    def toCnf(self) -> cnf.Cnf:
        return [[Literal(self.name())]]

class Negation(Formula):
    def __init__(self, orig : Formula) -> None:
        Formula.__init__(self, [orig])
    def originalFormula(self) -> Formula:
        return self.subf()[0]
    def isSatisfied(self, v : Valuation) -> bool:
        return not self.originalFormula().isSatisfied(v)
    def toString(self) -> str:
        return "-%s" % (self.originalFormula().toString())
    def toCnf(self) -> cnf.Cnf:
        forcnf = originalFormula(self).toCnf()
        if (forcnf.length() < 1):
            ret = forcnf;
        if (forcnf.length() == 1):
            retf = [];
            cla = forcnf[0];
            for lit in cla:
                retf.append(-lit)              
            ret = [retf]
        else: 
            retf = [];
            for lit1 in forcnf[0]:
                for lit2 in forcnf[1]:
                    retf.append( [(-lit1), (-lit2)])              
            ret = (Disjunction([retf]) + [forcnf[2:]]).toCnf()

        return ret;


class Disjunction(Formula):
    def __init__(self, subs : List[Formula]) -> None:
        Formula.__init__(self, subs)
    def isSatisfied(self, v : Valuation) -> bool:
        return any(f.isSatisfied(v) for f in self.subf())
    def toString(self) -> str:
        return '(' + '|'.join(f.toString() for f in self.subf()) + ')'
    def toCnf(self) -> cnf.Cnf:
        if (self.subf().length() < 2):
          ret = self.subf()
        else: 
          retf = []
          for f1 in self.items()[0]:
              for f2 in self.items()[1]:  
                  retf.append(toCnf([f1,f2]))              
          ret = retf+self.subf()[2:]
        return ret;

class Conjunction(Formula):
    def __init__(self, subs : List[Formula]) -> None:
        Formula.__init__(self, subs)
    def isSatisfied(self, v : Valuation) -> bool:
        return all(f.isSatisfied(v) for f in self.subf())
    def toString(self) -> str:
        return '(' + '&'.join(f.toString() for f in self.subf()) + ')'
    def toCnf(self) -> cnf.Cnf:
        v=[]
        for cla in self.subf():
            v.append(toCnf(cla))
        return v

class BinaryFormula(Formula):
    connective = ''
    def __init__(self, left : Formula, right : Formula) -> None:
        Formula.__init__(self, [left, right])
    def leftSide(self) -> Formula:
        return self.subf()[0]
    def rightSide(self) -> Formula:
        return self.subf()[1]
    def toString(self) -> str:
        return '(%s%s%s)' % (self.leftSide().toString(), self.connective, self.rightSide().toString())

class Implication(BinaryFormula):
    connective = '->'
    def isSatisfied(self, v : Valuation) -> bool:
        return (not self.leftSide().isSatisfied(v)) or self.rightSide().isSatisfied(v)
    def toCnf(self) -> cnf.Cnf:
        return Disjunction([Negation((self.leftSide()).toCnf()), (self.rightSide()).toCnf())).toCnf()

class Equivalence(BinaryFormula):
    connective = '<->'
    def isSatisfied(self, v : Valuation) -> bool:
        return self.leftSide().isSatisfied(v) == self.rightSide().isSatisfied(v)
    def toCnf(self) -> cnf.Cnf:
        return Conjunction(Disjunction([Negation((self.leftSide()).toCnf()),(self.rightSide()).toCnf())),
                           Disjunction([Negation((self.rightSide()).toCnf()), (self.leftSide()).toCnf()))).toCnf()

# vim: set sw=4 ts=4 sts=4 et :
