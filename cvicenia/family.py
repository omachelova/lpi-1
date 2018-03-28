#!/usr/bin/env python3

import os
import sys
#sys.path[0:0] = [os.path.join(sys.path[0], '.')]
sys.path[0:0] = [os.path.join(sys.path[0], '../../examples/sat')]

import string
import subprocess

from typing import List, Mapping, Union, TextIO
Valuation = Mapping[str, bool]

""" Nastavte na True ak chcete vypisat detaily riesenia """
Debbug_On = False

""" Nastavte na True ak chcete aby pouzival solver na dokazanie nemoznosti nastat """
proveusingsolver = False

class DimacsWriter(object):
    """ A helper class that writes clauses to a DIMACS format file. """
    def __init__(self, filename, mode = 'w'):
        """ Create a new writer that writes to *filename*.
            You can use ``mode='a'`` to append to an existing file
            instead of rewriting it.
        """
        self.fn = filename
        self.f = open(filename, mode)

    def filename(self):
        """ Returns the filename that this writer writes to as a string."""
        return self.fn

    def writeLiteral(self, lit):
        """ Writes a single literal (positive or negative integer).
            Use finishClause to finis this clause (write a zero).
        """
        self.f.write('{} '.format(lit))

    def finishClause(self):
        """" Finishes current clause (writes a zero).
            Note that if no clause was started (through *writeLiteral*),
            it will create an empty clause, which is always false!
        """
        self.f.write(' 0\n')
        self.f.flush()

    def writeClause(self, clause):
        """ Writes a single clause.
            *clause* must be a list of literals (positive or negative integers).
        """
        for l in clause:
            self.writeLiteral(l)
        self.finishClause()

    def writeImpl(self, left, right):
        """ Writes an implication *left* => *right*. """
        self.writeClause([-left, right])

    def writeComment(self, comment):
        """ Writes a comment.
            Note that this does not work inside an unfinished clause!
        """
        for line in comment.split('\n'):
            self.f.write('c {}\n'.format(line))

    def closed(self):
        """ Returs True if the output file has been already closed. """
        return self.f.closed

    def close(self):
        """ Closes the output file. """
        self.f.close()


class SatSolver(object):
    """ A helper class that manages SAT solver invocation. """

    def __init__(self, solverPath = None):
        """ Creates a new SAT solver.
            Use *solverPath* to specify an optional location where to look
            for SAT solver binary (it will be looked for in a set of default
            locations).
        """

        self.paths = []
        if solverPath:
            self.paths.append(solverPath)

        if sys.platform.startswith('linux'):
            self.paths += [
#                    'minisat', 'MiniSat_v1.14_linux',
                    './minisat', './MiniSat_v1.14_linux',
                    '../tools/lin/minisat',
                    '../../tools/lin/minisat'
                ]
        elif sys.platform.startswith('win'):
            self.paths += [
                    'miniexe', 'MiniSat_v1.14.exe',
                    '../tools/win/miniexe',
                    '../../tools/win/miniexe',
                ]
        else:
            pass # empty solver paths will fall back to try 'minisat'

        # default fall for all
        self.paths.append('minisat')

    def getSolverPath(self):
        """ Returns the path to solver binary. """
        for fn in self.paths:
            try:
                subprocess.check_output([fn, '--help'], stderr = subprocess.STDOUT)
                #sys.stderr.write('using sat solver:  "%s"\n' % fn)
                return fn
            except OSError:
                pass
        raise IOError('Solver executable not found!')

    def solve(self, theory, output):
        """ Use SAT solver to solve a theory, which is either the name
            of a file (in DIMACS format) or an instance of DimacsWriter.
            Writes the SAT solvers output to a file named *output*.
            Returns a tuple (sat, solution), where sat is True or False
            and solution is a list of positive or negative integers
            (an empty list if sat is False).
        """
        if isinstance(theory, DimacsWriter):
            if not theory.closed():
                theory.close()
            theory = theory.filename()

        try:
            self.output = subprocess.check_output(
                    [self.getSolverPath(), theory, output],
                    stderr = subprocess.STDOUT,

                    )
        except subprocess.CalledProcessError:
            # minisat has weird return codes
            pass

        with open(output) as f:
            sat = f.readline()
            if sat.strip() == 'SAT':
                sol = f.readline()
                return (
                        True,
                        [int(x) for x in sol.split()][:-1]
                )
            else:
                return (False, [])

class Literal(object):
    """ Reprezentacia literalu (premenna alebo negovana premenna) v CNF formule. """
    def __init__(self, name : str) -> None:
        """ Vytvori novy, kladny (nenegovany) literal pre premennu name. """
        self.name = name
        self.neg = False

    @staticmethod
    def Not(name : str) -> 'Literal':
        """ Vytvory novy, negovany literal pre premennu name. """
        lit = Literal(name)
        lit.neg = True
        return lit

    @staticmethod
    def fromInt(i : int, varMap : Union['VariableMap', Mapping[int, str]]) -> 'Literal':
        """ Vytvori novy literal podla ciselneho kodu a mapovania. """
        if isinstance(varMap, VariableMap):
            lit = Literal(varMap.reverse()[abs(i)])
        else:
            lit = Literal(varMap[abs(i)])
        if i < 0:
            lit.neg = True
        return lit

    def __neg__(self) -> 'Literal':
        """ Vrati novy literal, ktory je negaciou tohoto.

        Toto je specialna metoda, ktoru python zavola, ked
        na nejaky objekt pouzijeme operator -, t.j.:
        l1 = Literal('a')
        l2 = -l1  # l2 je teraz negaciou l1
        """
        lit = Literal(self.name)
        lit.neg = not self.neg
        return lit

    def toString(self) -> str:
        """ Vrati textovu reprezentaciu tohoto literalu. """
        if self.neg:
            return "-" + self.name
        else:
            return self.name

    def __str__(self) -> str:
        return self.toString()

    def __repr__(self) -> str:
        return self.__class__.__name__ + '(' + self.toString() + ')'


    def isSatisfied(self, v : Valuation) -> bool:
        """ Urci, ci je tento literal splneny ohodnotenim v. """
        return bool(bool(self.neg) ^ bool(v[self.name]))

    def writeToFile(self, outFile : TextIO, varMap : 'VariableMap'):
        """ Zapise literal do suboru outFile s pouzitim mapovania premennych varMap. """
        if self.neg:
            outFile.write('-%d' % varMap[self.name])
        else:
            outFile.write('%d' % varMap[self.name])

class VariableMap(dict):
    """ Mapovanie mien premennych na cisla.

    Premennym vzdy priraduje suvisly usek cisel 1..n.
    """

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self._varNumMax = 0

    def __missing__(self, key):
        self._varNumMax += 1
        self[key] = self._varNumMax
        return self._varNumMax

    def reverse(self) -> Mapping[int, str]:
        """ Vrati reverzne mapovanie ako jednoduchy slovnik z cisel na mena premennych. """
        rev = {}
        for k,v in self.items():
            rev[v] = k
        return rev

    def extend(self, what : Union[str, Literal, 'Cnf', 'Clause']):
        """ Rozsiri tuto mapu o premenne v danej cnf / klauze / literale. """
        if isinstance(what, str):
            self[what]
        if isinstance(what, Literal):
            self[what.name]
        elif isinstance(what, list): # Cnf and Clause are list-s
            for x in what:
                self.extend(x)

    def writeToFile(self, outFile : TextIO, prefix : str = ''):
        """ Zapise mapu do suboru outFile.

        Jedna premenna na jeden riadok (t.j. nefunguje na premenne
        s koncom riadku v nazve).

        Volitelny retazec prefix sa prida pred kazdy riadok
        (napriklad 'c ' vyrobi komentare pre dimacs cnf format).
        """
        outFile.write('%s%d\n' % (prefix, self._varNumMax))
        rev = self.reverse()
        for i in range(self._varNumMax):
            outFile.write('%s%s\n' % (prefix, rev[i+1]))

    @staticmethod
    def readFromFile(inFile : TextIO, prefix : str = '') -> 'VariableMap':
        """ Nacita novu mapu zo suboru inFile a vrati ju.

        Ak je uvedeny volitelny retazec prefix, tento sa odstrani
        zo zaciatku kazdeho riadku (napriklad 'c ' ak je mapa ulozena
        ako komentar v dimacs cnf formate).
        """
        def removePrefix(s : str) -> str:
            """Remove a prefix from string if present."""
            return s[len(prefix):] if s.startswith(prefix) else s

        varMap = VariableMap([])
        n = int(removePrefix(inFile.readline()))
        for i in range(n):
            varMap[removePrefix(inFile.readline().rstrip('\n'))]
        return varMap



class Clause(list):
    """ Reprezentacia klauzy (pole literalov). """
    def __init__(self, *args, **kwargs):
        """ Vytvori novu klauzu obsahujucu argumenty konstruktora. """
        list.__init__(self, *args, **kwargs)
        for lit in self:
            if not isinstance(lit, Literal):
                raise TypeError('Clause can contain only Literal-s')

    def toString(self) -> str:
        """ Vrati textovu reprezentaciu tejto klauzy (vid zadanie). """
        return ' '.join(lit.toString() for lit in self)

    def __str__(self) -> str:
        return self.toString()

    def isSatisfied(self, v : Valuation) -> bool:
        """ Urci, ci je tato klauza splnena ohodnotenim v. """
        for lit in self:
            if lit.isSatisfied(v):
                return True
        return False

    def writeToFile(self, oFile : TextIO, varMap : Mapping[str, int]):
        """ Zapise klauzu do suboru outFile v DIMACS formate
            pricom pouzije varMap na zakodovanie premennych na cisla.

        Klauzu zapise na jeden riadok (ukonceny znakom konca riadku).
        """
        for lit in self:
            lit.writeToFile(oFile, varMap)
            oFile.write(' ')
        oFile.write(' 0\n')

    @staticmethod
    def readFromFile(inFile : TextIO, varMap : VariableMap) -> 'Clause':
        """ Nacita novu klauzu zo suboru inFile a vrati ju ako vysledok.
        Mozete predpokladat, ze klauza je samostatne na jednom riadku.
        Ak sa z aktualneho riadku na vstupe neda nacitat korektna klauza,
        vyhodi vynimku IOError.
        """
        line = inFile.readline()
        if line is None:
            raise IOError('End of file')

        rVarMap = varMap.reverse()

        cls = Clause()
        ints = [int(x) for x in line.split()]
        if len(ints) < 1 or ints[-1] != 0:
            raise IOError('Bad clause')
        for i in ints[:-1]:
            if i == 0:
                raise IOError('Bad clause (0 inside)')
            elif i < 0:
                cls.append(Literal.Not(rVarMap[abs(i)]))
            else:
                cls.append(Literal(rVarMap[abs(i)]))
        return cls

class Cnf(list):
    """ Reprezentacia Cnf formuly ako pola klauz. """
    def __init__(self, *args, **kwargs):
        """ Vytvori novu Cnf formulu obsahujucu argumenty konstruktora. """
        list.__init__(self, *args, **kwargs)
        for cls in self:
            if not isinstance(cls, Clause):
                raise TypeError('Cnf can contain only Clause-s')

    def toString(self) -> str:
        """ Vrati textovu reprezentaciu tejto formuly (vid zadanie). """
        return ''.join([ cls.toString() + '\n' for cls in self])

    def __str__(self) -> str:
        return self.toString()

    def isSatisfied(self, v : Valuation) -> bool:
        """ Urci, ci je tato formula splnena ohodnotenim v. """
        for cls in self:
            if not cls.isSatisfied(v):
                return False
        return True

    def extendVarMap(self, varMap : VariableMap):
        """ Rozsiri varMap o premenne v tejto formule. """
        for cls in self:
            cls.extendVarMap(varMap)

    def writeToFile(self, oFile : TextIO, varMap : VariableMap):
        """ Zapise klauzu do suboru outFile v DIMACS formate
            pricom pouzije varMap na zakodovanie premennych na cisla
            a zapise kazdu klauzu na jeden riadok.
        """
        for cls in self:
            cls.writeToFile(oFile, varMap)

    @staticmethod
    def readFromFile(inFile : TextIO, varMap : VariableMap) -> 'Cnf':
        """ Nacita novu formulu zo suboru inFile a vrati ju ako vysledok.
        Mozete predpokladat, ze kazda klauza je samostatne na jednom riadku.
        """
        cnf = Cnf()
        while True:
            try:
                cls = Clause.readFromFile(inFile, varMap)
            except IOError:
                break
            append(cls)
        return cnf


Valuation = Mapping[str, bool]

#Literal
# Clause
#  Cnf
#   Formula
#    BinaryFormula
#     Implication
#     Equivalence

class Clause2(Clause):
    def __init__(self, clauseorig : Clause) -> None:
        Clause.__init__(self, clauseorig)
        self._clauseorig = Clause(clauseorig)
    def spakuj(self, vratprazdno:bool) -> Clause:
        retcla = Clause()
        for i in range(len(self._clauseorig)):
            name_i = self._clauseorig[i].name
            if (name_i[0:1] == "-"):
                neg_i = 1
                name_i = name_i[1:]
            else:
                neg_i = 0

            j = 0 
            while (j < len(retcla)):
                name_j = retcla[j].name
                if (name_j[0:1] == "-"):
                    neg_j = 1
                    name_j = name_j[1:]
                else:
                    neg_j = 0

                if ((name_i < name_j) or ((name_i == name_j) and (neg_i <= neg_j))):
                    break
                j+=1
            if (j < len(retcla)):
                if (name_i == name_j):
                    if (neg_i == neg_j):
                        continue
            retcla2 = retcla[0:j] + [self._clauseorig[i]] + retcla[j:]
            retcla = retcla2               

        name_new = ""
        neg_new = 2

        for i in range(len(retcla)):
            name_old = name_new
            neg_old = neg_new

            name_new = retcla[i].name
            if (name_new[0:1] == "-"):
                neg_new = 1
                name_new = name_new[1:]
            else:
                neg_new = 0

            if ((name_new == name_old) and (neg_new != neg_old)):
                if vratprazdno:
                    retcla = Clause()
                    return retcla
  
        return retcla

class Formula(object):
    def __init__(self, subs : List['Formula']= []) -> None:
        self.m_subf = subs # type: List[Formula]
    def subf(self) -> List['Formula']:
        return self.m_subf
    def isSatisfied(self, v : Mapping[str,bool]) -> bool:
        return False
    def toString(self) -> str:
        return "INVALID"
    def convertfromlist(self, listret) -> Cnf:
        return listret
    def lenilja(self) -> int:
        return 0
    def __str__(self) -> str:
        return self.toString()
    def __repr__(self) -> str:
        return self.__class__.__name__ + '(' + ','.join([ repr(f) for f in self.subf()]) + ')'

    def toCnfInternal(self) -> Cnf:
        return Cnf()

    def toCnf(self) -> Cnf:
        ftheory = Cnf(self.toCnfInternal())
        ftheorycnf = self.provecnf(ftheory)

        if (len(ftheory) == 0):
            retf = Cnf()

            lit1 = Literal("a")
            lit1neg = Literal("a")
            lit1neg = -lit1neg
            clause1neg =    Clause([lit1, lit1neg])
            retf1= Cnf([clause1neg])
            return retf1

        if (self.provealwaysfalse(ftheory, ftheorycnf)):
            retf = Cnf()

            lit1 = Literal(ftheory[0][0].name)
            lit1neg = Literal(ftheory[0][0].name)
            lit1neg = -lit1neg
            clause1 =    Clause([lit1])
            clause1neg = Clause([lit1neg])
            retf1= Cnf([clause1, clause1neg])
            return retf1
        else:
            return ftheorycnf

    def provecnf(self, teoria:Cnf) -> Cnf:
        c = teoria

        c2 = Cnf()
        for cla in c:
            cla2 = Clause()

            for lit in cla:
                if (lit.name[0:1] == "-"):
                    lit2 = Literal(lit.name[1:])
                    lit2 = -lit2
                else:
                    lit2 = Literal(lit.name)
                cla2.append(lit2)
            c2.append(cla2)

        return c2

    def provealwaysfalse(self, teoriavzdynepravda:Cnf, teoria:Cnf) -> bool:
        if not proveusingsolver:
            return False

        c = teoriavzdynepravda
        c2= teoria

        varMap = VariableMap()
        varMap.extend(c2)

        solver = SatSolver()
        w = DimacsWriter("provealwaysfalse_problem_cnf_in.txt")

        for clapom in c:

            clause=Clause2(clapom)
            cla = Clause2.spakuj(clause, True)
            if len(cla) == 0:
                continue

            for lit in cla:
                if (lit.name[0:1] == "-"):
                    w.writeLiteral(varMap[lit.name[1:]] + len(varMap))
                else:
                    w.writeLiteral(varMap[lit.name])
            w.finishClause()

        i = 0
        while (i < len(varMap)):
            i+=1
            w.writeLiteral((i + len(varMap)))
            w.writeLiteral(i)
            w.finishClause()

            w.writeLiteral(-(i + len(varMap)))
            w.writeLiteral(-(i))
            w.finishClause()

        w.close()
        ok, sol = solver.solve(w, "provealwaysfalse_problem_cnf_out.txt")

        if ok:
            #for x in sol:
            #    if x>0:
            #print("NIE VZDY LOZ!!!!!!!")
            return False
        #print("VZDY LOZ!!!!!!!")
        return True

class Variable(Formula, Cnf):
    def __init__(self, name : str) -> None:
        Formula.__init__(self)
        self._name = name
    def name(self) -> str:
        return self._name
    def isSatisfied(self, v : Mapping[str, bool]) -> bool:
        return v[self.name()]
    def toString(self) -> str:
        return self.name()
    def lenilja(self) -> int:
        return 1
    def __repr__(self) -> str:
        return "Variable(%r)" % (self.name(),)
    def toCnfInternal(self) -> Cnf:
        retf1 = Cnf([Clause([Literal(self.name())])])
        retf2 = Cnf(retf1)
        return retf2

class Negation(Formula, Cnf):
    def __init__(self, orig : Formula) -> None:
        Formula.__init__(self, [orig])
    def originalFormula(self) -> Formula:
        return self.subf()[0]
    def isSatisfied(self, v : Valuation) -> bool:
        return not self.originalFormula().isSatisfied(v)
    def toString(self) -> str:
        return "-%s" % (self.originalFormula().toString())
    def toCnfInternal(self) -> Cnf:
        forcnf = (self.subf()[0]).toCnfInternal()
        if (len(forcnf) == 0):
            retf = Cnf()
            lit    = Literal("A")
            litneg = Literal("-A")

            clause1 = Clause()
            clause1.append(lit)

            clause2 = Clause()
            clause2.append(litneg)

            retf.append(clause1)
            retf.append(clause2)

            retf1 = Cnf(retf)
            ret = self.convertfromlist(retf1)
        elif (len(forcnf) >= 1):
            retf = [""]
            k = 0
            while (k < len(forcnf)):
                retfold = retf
                retf = Cnf()
                for cla in retfold:
                    for lit1 in forcnf[k]:
                        if (lit1.name[0:1] == "-"):
                            lit2 = Literal(lit1.name[1:])
                        else:
                            lit2 = Literal("-"+lit1.name)
                        if k == 0:
                            clause2 = Clause([lit2])
                        else:
                            clause2 = Clause(cla)
                            clause2.append(lit2)
                        lit12neg = Clause(clause2)
                        clause=Clause2(lit12neg)
                        retfs = Clause2.spakuj(clause, True)
                        if len(retfs) > 0:
                            retf.append(Clause(retfs))
                    if (k == 0):
                        break
                k=k+1
            retf1 = Cnf(retf)
            ret = self.convertfromlist(retf1)
        return ret

class Disjunction(Formula, Cnf):
    def __init__(self, subs : List[Formula]) -> None:
        Formula.__init__(self, subs)
    def isSatisfied(self, v : Valuation) -> bool:
        return any(f.isSatisfied(v) for f in self.subf())
    def toString(self) -> str:
        return '(' + '|'.join(f.toString() for f in self.subf()) + ')'
    def toCnfInternal(self) -> Cnf:

        if (len(self.subf()) < 1):
            lit1 = Literal("a")
            lit1neg = Literal("-a")
            clause1    =    Clause([lit1])
            clause1neg =    Clause([lit1neg])
            ret= Cnf([clause1, clause1neg])
        elif (len(self.subf()) == 1):
            forcnf = (self.subf()[0]).toCnfInternal()
            ret = forcnf
        else:
            forcnf = (self.subf()[0]).toCnfInternal()
            forcnf1 =    Cnf(
                                 (Disjunction(self.subf()[1:])
                                 ).toCnfInternal()
                                )
            retf = Cnf()
            for cla1 in forcnf:
                for cla2 in forcnf1:
                    cla12= Clause(cla1 + cla2);
                    clause=Clause2(cla12);
                    retfs = Clause2.spakuj(clause, True)
                    if len(retfs) > 0:
                        retf.append(Clause(retfs))
            retf1 = Cnf(retf)
            ret = self.convertfromlist(retf1)

        return ret

class Conjunction(Formula, Cnf):
    def __init__(self, subs : List[Formula]) -> None:
        Formula.__init__(self, subs)
    def isSatisfied(self, v : Valuation) -> bool:
        return all(f.isSatisfied(v) for f in self.subf())
    def toString(self) -> str:
        return '(' + '&'.join(f.toString() for f in self.subf()) + ')'
    def toCnfInternal(self) -> Cnf:
        if (len(self.subf()) < 1):
            lit1 = Literal("a")
            lit1neg = Literal("-a")
            clause1neg =    Clause([lit1, lit1neg])
            ret= Cnf([clause1neg])
            return ret
        retf = Cnf()
        for f1 in self.subf():
            f2 = f1.toCnfInternal()
            for lit1 in f2:
                clause=Clause2(lit1)
                retfs = Clause2.spakuj(clause, True)
                if len(retfs) > 0:
                    retf.append(Clause(retfs))
        retf1 = (Cnf(retf))
        retf2 = self.convertfromlist(retf1)
        return retf2

class BinaryFormula(Formula, Cnf):
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
    def toCnfInternal(self) -> Cnf:
        return                   (Disjunction(
                                              [
                                               Negation((self.leftSide()))
                                               , 
                                               (self.rightSide())
                                              ]
                                             )).toCnfInternal()

class Equivalence(BinaryFormula):
    connective = '<->'
    def isSatisfied(self, v : Valuation) -> bool:
        return self.leftSide().isSatisfied(v) == self.rightSide().isSatisfied(v)
    def toCnfInternal(self) -> Cnf:
        return (
                 Conjunction(
                                 [
                                  Disjunction(
                                              [
                                               Negation((self.leftSide()))
                                               , 
                                               (self.rightSide())
                                              ]
                                             )
                                  ,
                                  Disjunction(
                                              [Negation((self.rightSide()))
                                               , 
                                               (self.leftSide())
                                              ]
                                             ) 
                                ]
                            )
                  ).toCnfInternal()

def find_is_possible(theory):
    """ Zisti ci theory je mozne. """

    f = theory
        
    c = f.toCnf()

    varMap = VariableMap()
    varMap.extend(c)

    fn = "problem_possible.txt"
    of = open(fn, "w")
    c.writeToFile(of, varMap)
    of.close()

    solver = SatSolver()
    ok, sol = solver.solve(fn, "problem_possible_out.txt")
    if ok:
        if Debbug_On:
            print("!!! {}   JE TO MOZNE   {} lebo:".format(
                theory.toString(),
                ', '.join([x.toString() for x in theory])))
            print("  {}".format(repr(sol)))
            revMap = varMap.reverse()
            print("  {}".format(repr([str(Literal.fromInt(i,varMap)) for i in sol])))

        """
        Virginia, Dorothy, George, Howard
        V, D, G, H
        1, 2, 3, 4
        """

        for i in sol:
            if (repr(str(Literal.fromInt(i,varMap))) == "'j14'"):
                print("otec: Virginia")
            elif (repr(str(Literal.fromInt(i,varMap))) == "'j24'"):
                print("otec: Dorothy")
            elif (repr(str(Literal.fromInt(i,varMap))) == "'j34'"):
                print("otec: George")
            elif (repr(str(Literal.fromInt(i,varMap))) == "'j44'"):
                print("otec: Howard")

        for i in sol:
            if (repr(str(Literal.fromInt(i,varMap))) == "'j13'"):
                print("matka: Virginia")
            elif (repr(str(Literal.fromInt(i,varMap))) == "'j23'"):
                print("matka: Dorothy")
            elif (repr(str(Literal.fromInt(i,varMap))) == "'j33'"):
                print("matka: George")
            elif (repr(str(Literal.fromInt(i,varMap))) == "'j43'"):
                print("matka: Howard")

        for i in sol:
            if (repr(str(Literal.fromInt(i,varMap))) == "'j12'"):
                print("syn: Virginia")
            elif (repr(str(Literal.fromInt(i,varMap))) == "'j22'"):
                print("syn: Dorothy")
            elif (repr(str(Literal.fromInt(i,varMap))) == "'j32'"):
                print("syn: George")
            elif (repr(str(Literal.fromInt(i,varMap))) == "'j42'"):
                print("syn: Howard")

        for i in sol:
            if (repr(str(Literal.fromInt(i,varMap))) == "'j11'"):
                print("dcera: Virginia")
            elif (repr(str(Literal.fromInt(i,varMap))) == "'j21'"):
                print("dcera: Dorothy")
            elif (repr(str(Literal.fromInt(i,varMap))) == "'j31'"):
                print("dcera: George")
            elif (repr(str(Literal.fromInt(i,varMap))) == "'j41'"):
                print("dcera: Howard")

    else:
        if Debbug_On:
            print("{}   NIE JE TO MOZNE   {}".format(
                theory.toString(),
                ', '.join([x.toString() for x in theory])))

        print("NIE JE TO MOZNE")


Not = Negation
Var = Variable
And = Conjunction
Or = Disjunction
Impl = Implication
Equ = Equivalence


#find_is_possible(
#  And([Var('v11'), Not(Var('v11'))])
#)
#find_is_possible(
#  And([Var('v11'), Not(Var('v12'))])
#)

"""
s{Osoba} JE STARSIA ASPON AKO {Osoba}
{s11, s12, s13, s14, s21, s22, s23, s24, s31, s32, s33, s34, s41, s42, s43, s44}
{Osoba}
Virginia, Dorothy, George, Howard
V, D, G, H
1, 2, 3, 4

j{Osoba}JE{Rodinny prislusnik}
{j11, j12, j13, j14, j21, j22, j23, j24, j31, j32, j33, j34, j41, j42, j43, j44}

{Rodinny prislusnik}
dcera, syn, matka, otec
d, s, m, o
1, 2, 3, 4
"""


"""
alebo jedno je vacsie rovne ako druhe, alebo druhe je vacsie rovne ako prve
&& {sij || sji}

TRANZITIVNOST NEROVNOSTI VEKOV (NEOSTREJ)

sij && sjk => sik

&& {!sij || !sjk || sik}
"""

konzistencia_usporiadania_veku_zoz=[]

for i in range(1, 5):
    for j in range(1, 5):
        konzistencia_usporiadania_veku_zoz.append( Or( [Var('s'+ ('%d' % (i)) + ('%d' % (j)) ), 
                                                        Var('s'+ ('%d' % (j)) + ('%d' % (i)) )]
                                                      ) )

for i in range(1, 5):
    for j in range(1, 5):
        for k in range(1, 5):
            konzistencia_usporiadania_veku_zoz.append( Or( [
                                                            Not(Var('s'+ ('%d' % (i)) + ('%d' % (j)) )),
                                                            Not(Var('s'+ ('%d' % (j)) + ('%d' % (k)) )),
                                                            Var('s'+ ('%d' % (i)) + ('%d' % (k)) )
                                                           ]
                                                          )
                                                      )
                                                       
konzistencia_usporiadania_veku = And(konzistencia_usporiadania_veku_zoz)

konzistencia_rozmiestnenia_rodiny_zoz=[]

for i in range(1, 5):
    konzistencia_rozmiestnenia_rodiny_zoz.append( Or(
                                                     [ 
                                                      Var('j'+ ('%d' % (i)) + ('%d' % (1)) ), 
                                                      Var('j'+ ('%d' % (i)) + ('%d' % (2)) ), 
                                                      Var('j'+ ('%d' % (i)) + ('%d' % (3)) ), 
                                                      Var('j'+ ('%d' % (i)) + ('%d' % (4)) )
                                                     ]
                                                    )
                                                 ) 
    konzistencia_rozmiestnenia_rodiny_zoz.append( Or( 
                                                     [
                                                      Var('j'+ ('%d' % (1)) + ('%d' % (i)) ), 
                                                      Var('j'+ ('%d' % (2)) + ('%d' % (i)) ), 
                                                      Var('j'+ ('%d' % (3)) + ('%d' % (i)) ), 
                                                      Var('j'+ ('%d' % (4)) + ('%d' % (i)) )
                                                     ]
                                                    )
                                                 ) 

for i in range(1, 5):
    for j in range(1, 5):
        for k in range(1, 5):
            if (j != k):
                konzistencia_rozmiestnenia_rodiny_zoz.append( Or( 
                                                                 [
                                                                  Not(  Var('j'+ ('%d' % (i)) + ('%d' % (j)) ) ), 
                                                                  Not(  Var('j'+ ('%d' % (i)) + ('%d' % (k)) ) ) 
                                                                 ]
                                                                )
                                                             )

for i in range(1, 5):
    for j in range(1, 5):
        for k in range(1, 5):
            if (j != k):
                konzistencia_rozmiestnenia_rodiny_zoz.append( Or( 
                                                                 [
                                                                  Not(  Var('j'+ ('%d' % (j)) + ('%d' % (i)) ) ), 
                                                                  Not(  Var('j'+ ('%d' % (k)) + ('%d' % (i)) ) ) 
                                                                 ]
                                                                )
                                                             )

konzistencia_rozmiestnenia_rodiny = And(konzistencia_rozmiestnenia_rodiny_zoz)

"""
j{Osoba}JE{Rodinny prislusnik}
{j11, j12, j13, j14, j21, j22, j23, j24, j31, j32, j33, j34, j41, j42, j43, j44}

{Osoba}
Virginia, Dorothy, George, Howard
V, D, G, H
1, 2, 3, 4

{Rodinny prislusnik}
dcera, syn, matka, otec
d, s, m, o
1, 2, 3, 4

vylucenie zlych pohlavnych prislusnosti v klasickej Americkej rodine

{!Vs && !Vo}
{!Ds && !Do}
{!Gd && !Gm}
{!Hd && !Hm}

"""

konzistencia_pohlavi_rodiny_zoz=[]

for i in range(1, 5):
    for j in range(1, 5):
        if ( ((((int)((i-1) / 2))) % 2) != ((j-1) % 2) ):
            konzistencia_pohlavi_rodiny_zoz.append( Not(  Var('j'+ ('%d' % (i)) + ('%d' % (j)) ) )
                                                             )
konzistencia_pohlavi_rodiny = And(konzistencia_pohlavi_rodiny_zoz)

"""

vylucenie vekovych anomalii v klasickej Americkej rodine 

dieta nemoze byt starsie alebo rovne vekovo ako rodic

((Osoba{j1} == {syn alebo dcera}) && (Osoba{j2} == {mama alebo otec})) => !s{j1}{j2}

(((Osoba{j1} == dcera) alebo (Osoba{j1} == syn)) && ( (Osoba{j2} == mama) alebo (Osoba{j2} == otec})) => !s{j1}{j2}


(((Osoba{j1} == dcera) alebo (Osoba{j1} == syn)) && ( (Osoba{j2} == mama) alebo (Osoba{j2} == otec})) => !s{j1}{j2}


( (j{j1}1 || j{j1}2) && (j{j2}3 || j{j2}4) ) => !s{j1}{j2}



( (!j{j1}1 && !j{j1}2) || (!j{j2}3 && !j{j2}4) ) || !s{j1}{j2}


 (!j{j1}1 || !j{j2}3 || !s{j1}{j2}) && 
 (!j{j1}1 || !j{j2}4 || !s{j1}{j2}) && 
 (!j{j1}2 || !j{j2}3 || !s{j1}{j2}) && 
 (!j{j1}2 || !j{j2}4 || !s{j1}{j2}) 


"""

konzistencia_vekov_rodiny_zoz=[]

for j1 in range(1, 5):
    for j2 in range(1, 5):
        konzistencia_vekov_rodiny_zoz.append( Or([
                                                  Not(  Var('j'+ ('%d' % (j1)) + ('%d' % (1 )) ) ),
                                                  Not(  Var('j'+ ('%d' % (j2)) + ('%d' % (3 )) ) ),
                                                  Not(  Var('s'+ ('%d' % (j1)) + ('%d' % (j2)) ) )
                                                 ])
                                             )
        konzistencia_vekov_rodiny_zoz.append( Or([
                                                  Not(  Var('j'+ ('%d' % (j1)) + ('%d' % (1 )) ) ),
                                                  Not(  Var('j'+ ('%d' % (j2)) + ('%d' % (4 )) ) ),
                                                  Not(  Var('s'+ ('%d' % (j1)) + ('%d' % (j2)) ) )
                                                 ])
                                             )
        konzistencia_vekov_rodiny_zoz.append( Or([
                                                  Not(  Var('j'+ ('%d' % (j1)) + ('%d' % (2 )) ) ),
                                                  Not(  Var('j'+ ('%d' % (j2)) + ('%d' % (3 )) ) ),
                                                  Not(  Var('s'+ ('%d' % (j1)) + ('%d' % (j2)) ) )
                                                 ])
                                             )
        konzistencia_vekov_rodiny_zoz.append( Or([
                                                  Not(  Var('j'+ ('%d' % (j1)) + ('%d' % (2 )) ) ),
                                                  Not(  Var('j'+ ('%d' % (j2)) + ('%d' % (4 )) ) ),
                                                  Not(  Var('s'+ ('%d' % (j1)) + ('%d' % (j2)) ) )
                                                 ])
                                             )

konzistencia_vekov_rodiny = And(konzistencia_vekov_rodiny_zoz)

"""
osoby i1 a i2 su pokrvni pribuzni

!(i1 a i2 NIE SU pokrvni pribuzni)

to, ze 2 ludia NIE SU pokrvni pribuzni v klasickej americkej rodine je mozne prave vtedy, ak su to otec a matka
teda
!( (j{i1}3 && j{i2}4) || (j{i1}4 && j{i2}3) )   

( (!j{i1}3 || !j{i2}4) && (!j{i1}4 || !j{i2}3) )   

E1 <=> osoby "2" a "3" su pokrvni pribuzni

E1 <=> ( (!j23 || !j34) && (!j24 || !j33) )   

"""


definicie_E1_E4_zoz = []

definicie_E1_E4_zoz.append( 
                            Equ(
                                Var('E1'),  
                                And([
                                     Or([Not(Var('j23')), Not(Var('j34'))]),
                                     Or([Not(Var('j24')), Not(Var('j33'))])
                                    ])
                                )
                           )
 
definicie_E1_E4_zoz.append( 
                            Equ(
                                Var('E2'),  
                                Var('s43')
                               )
                           )

definicie_E1_E4_zoz.append( 
                            Equ(
                                Var('E3'),  
                                Not(Var('s14'))
                               )
                           )

definicie_E1_E4_zoz.append( 
                            Equ(
                                Var('E4'),  
                                Var('s12')
                               )
                           )

definicie_E1_E4 = And(definicie_E1_E4_zoz)


"""
prave 2 z E1, E2, E3, E4 su pravdive

{{!E1 || !E2 || !E3} && {!E1 || !E2 || !E4} && {!E1 || !E3 || !E4} && {!E2 || !E3 || !E4}}
===========================================================================================

{{E1 || E2 || E3} && {E1 || E2 || E4} && {E1 || E3 || E4} && {E2 || E3 || E4}}

!{{E1 && E2 && E3} || {E1 && E2 && E4} || {E1 && E3 && E4} || {E2 && E3 && E4}}


!{{!E1 && !E2 && !E3} || {!E1 && !E2 && !E4} || {!E1 && !E3 && !E4} || {!E2 && !E3 && !E4}}

{{E1 || E2 || E3} && {E1 || E2 || E4} && {E1 || E3 || E4} && {E2 || E3 || E4}}
===============================================================================




!{E1 && E2 && E3} || 

NEMAME 

===============================================================================
{{!E1 || !E2 || !E3} && {!E1 || !E2 || !E4} && {!E1 || !E3 || !E4} && {!E2 || !E3 || !E4}}
&&
{{E1 || E2 || E3} && {E1 || E2 || E4} && {E1 || E3 || E4} && {E2 || E3 || E4}}
===============================================================================

"""

pocet_pravdivych_E1_E4_zoz = []

pocet_pravdivych_E1_E4_zoz.append( 
                                  Or([Not(Var('E1')), Not(Var('E2')), Not(Var('E3'))
                                     ])
                                 )
pocet_pravdivych_E1_E4_zoz.append( 
                                  Or([Not(Var('E1')), Not(Var('E2')), Not(Var('E4'))
                                     ])
                                 )
pocet_pravdivych_E1_E4_zoz.append( 
                                  Or([Not(Var('E1')), Not(Var('E3')), Not(Var('E4'))
                                     ])
                                 )
pocet_pravdivych_E1_E4_zoz.append( 
                                  Or([Not(Var('E2')), Not(Var('E3')), Not(Var('E4'))
                                     ])
                                 )
 

pocet_pravdivych_E1_E4_zoz.append( 
                                  Or([Var('E1'), Var('E2'), Var('E3')
                                     ])
                                 )
pocet_pravdivych_E1_E4_zoz.append( 
                                  Or([Var('E1'), Var('E2'), Var('E4')
                                     ])
                                 )
pocet_pravdivych_E1_E4_zoz.append( 
                                  Or([Var('E1'), Var('E3'), Var('E4')
                                     ])
                                 )
pocet_pravdivych_E1_E4_zoz.append( 
                                  Or([Var('E2'), Var('E3'), Var('E4')
                                     ])
                                 )

pocet_pravdivych_E1_E4 = And(pocet_pravdivych_E1_E4_zoz)


if Debbug_On:
    find_is_possible(
      And(
          [
           konzistencia_usporiadania_veku,
           konzistencia_rozmiestnenia_rodiny,
           konzistencia_pohlavi_rodiny,
           konzistencia_vekov_rodiny,
           definicie_E1_E4,
           pocet_pravdivych_E1_E4,
           Not(And([Var('j13'), Var('j21'), Var('j42'),Var('j34')]))
          ]      
        )
    )

find_is_possible(
  And(
      [
       konzistencia_usporiadania_veku,
       konzistencia_rozmiestnenia_rodiny,
       konzistencia_pohlavi_rodiny,
       konzistencia_vekov_rodiny,
       definicie_E1_E4,
       pocet_pravdivych_E1_E4
      ]      
    )
)
  

# vim: set sw=4 ts=4 sts=4 et :
