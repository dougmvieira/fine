from collections import namedtuple
from functools import reduce
from itertools import chain
from sympy import Symbol, Function, Eq
from sympy.codegen.ast import CodeBlock, Assignment, Variable, complex128


FunctionTree = namedtuple("FunctionTree", ["symb", "func", "args", "expr"])

def tree_constr(symb_str, args, expr, func=None):
    symb = Symbol(symb_str)
    func = func if func else Function(symb_str)(
            *(arg.symb if type(arg) is FunctionTree else arg for arg in args))
    return FunctionTree(symb, func, args, expr)

def eq_form(tree):
    return Eq(tree.func, tree.expr, evaluate=False)

def modify_expr(func):
    def closure(tree):
        return FunctionTree(tree.symb, tree.func, tree.args, func(tree.expr))
    return closure

def diff(tree, arg):
    arg_symb = arg.symb if hasattr(arg, "symb") else arg
    diff_symb = Symbol("{}_{}".format(str(tree.symb), str(arg_symb)))
    diff_func = tree.func.diff(arg_symb)
    diff_expr = tree.expr.diff(arg_symb)
    return FunctionTree(diff_symb, diff_func, tree.args, diff_expr)

def chain_preserve_uniqueness_step(acc, it):
    return list(chain(acc, filter(lambda el: el not in acc, it)))

def chain_preserve_uniqueness(*its):
    return reduce(chain_preserve_uniqueness_step, its, [])

def traverse_expr(expr, func):
    args = [traverse_expr(arg, func) if type(arg) is not Symbol
            else arg for arg in expr.args]
    return func(expr, args)

def args_in_expr_func(_, args):
    acc_args = filter(lambda arg: type(arg) is list, args)
    just_args = filter(lambda arg: type(arg) is not list, args)
    return chain_preserve_uniqueness(just_args, *acc_args)

def args_in_expr(expr):
    return ([expr] if type(expr) is Symbol
            else traverse_expr(expr, args_in_expr_func))

def traverse_tree(tree, func):
    args = [traverse_tree(arg, func) if type(arg) is FunctionTree
            else arg for arg in tree.args]
    return func(tree, args)

def prepend_precedents_func(tree, args):
    precedents = chain_preserve_uniqueness(
        *filter(lambda arg: type(arg) is list, args))
    return chain_preserve_uniqueness(precedents, [tree])

def prepend_precedents(tree):
    return traverse_tree(tree, prepend_precedents_func)

def diffs_func(tree, args):
    node_diffs = (diff(tree, arg) for arg in tree.func.args)
    args_diffs = chain_preserve_uniqueness(
        *filter(lambda arg: type(arg) is list, args))
    return chain_preserve_uniqueness(args_diffs, node_diffs)

def diffs(tree):
    return traverse_tree(tree, diffs_func)

def chain_rule_func(symb):
    def closure(tree, args):
        dxs = [arg.diff(symb) if type(arg) is Symbol else arg for arg in args]
        f_primes = [tree.func.diff(arg) for arg in tree.func.args]
        return sum(f_prime*dx for f_prime, dx in zip(f_primes, dxs))
    return closure

def chain_rule(tree, symb):
    return traverse_tree(tree, chain_rule_func(symb))

def chain_rule_tree(tree, tree_diffs, *symbs):
    ret_symbs = (Symbol("{}_{}_chain".format(str(tree.symb), str(symb)))
                 for symb in symbs)
    funcs = (tree.func.func(*symbs).diff(symb) for symb in symbs)
    exprs = [func_to_symbs(chain_rule(tree, symb), tree_diffs)
             for symb in symbs]
    argss = ([t for t in tree_diffs if t.symb in args_in_expr(expr)]
             for expr in exprs)
    return [FunctionTree(symb, func, args, expr)
            for symb, func, args, expr in zip(ret_symbs, funcs, argss, exprs)]

def func_to_symbs(expr, trees):
    return expr.subs((tree.func, tree.symb) for tree in trees)

def tree_codeblock(*trees):
    tree_sequence = chain_preserve_uniqueness(*map(prepend_precedents, trees))
    assigments = (Assignment(tree.symb, tree.expr) for tree in tree_sequence)
    declarations = (Variable(tree.symb, type=complex128).as_Declaration()
                    for tree in tree_sequence)
    return CodeBlock(*declarations, *assigments)
