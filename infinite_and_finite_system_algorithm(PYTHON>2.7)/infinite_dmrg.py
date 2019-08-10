# Simple DMRG code.
# This code is a basic implementation of the infinite system algorithm.
#
# Author: William Gabriel Carreras Oropesa.
# Institution: Atomic Center of Bariloche.


# This code is compatible with any python version greater than 2.6,
# The following line provides consistency between python_2 and python_3.
from __future__ import print_function, division

import numpy as np
from scipy.sparse import kron, identity
from scipy.sparse.linalg import eigsh

from collections import namedtuple


Block = namedtuple("Block", ["length", "basis_size", "operator_dict"])
Enlarged_Block = namedtuple("Enlarged_Block",  ["length", "basis_size", "operator_dict"])


def is_valid_block(block):
    for op in block.operator_dict.values():
        if op.shape[0] != block.basis_size or op.shape[1] != block.basis_size:
            return False
        else:
            return True


is_valid_enlarged_block = is_valid_block

# Model specific code for the Heisenberg XXZ chain.
model1_d = 2    # single-site basis size

Sz1 = np.array([[0.5, 0], [0, -0.5]], dtype="d")    # single-site S^z operator.
Sp1 = np.array([[0, 1], [0, 0]], dtype="d")         # single-site S^+ operator.

H1 = np.array([[0, 0], [0, 0]], dtype="d")          # single-site, H is zero.


def H2(Sz1, Sp1, Sz2, Sp2):     # two-sites part of H.
    """
    This function take the operators of one sites, of tow different Hilbert spaces,
    and creates the operators of tow sites.
    :param Sz1: projection of spin in Oz, site 1 (np.array).
    :param Sp1: operator of entanglement in site 1 (np.array).
    :param Sz2: projection of spin in Oz, site 2 (np.array).
    :param Sp2: operator of entanglement in site 2 (np.array).
    :return:  Hamiltonian operator of tow sites.   (np.array).
    """
    part_1 = 0.5 * kron(Sp1, Sp2.conjugate().transpose())
    part_2 = 0.5 * kron(Sp1.conjugate().transpose(), Sp2)
    part_3 = kron(Sz1, Sz2)
    return part_1 + part_2 + part_3


# Here we defines the initial block.
initial_block = Block(length=1, basis_size=model1_d, operator_dict={
    "H": H1,
    "conn_Sz": Sz1,
    "conn_Sp": Sp1,
})


def enlarge_block(block):
    """
    This function enlarges the provided block in a single site.
    :param block: block to enlarge (dictionary)
    :return: Enlarger block (dictionary)
    """
    my_block = block.basis_size
    o = block.operator_dict

    # Create the new operators for the enlarger block
    enlarge_block_dict = {
        "H": kron(o["H"], identity(model1_d)) + kron(identity(my_block), H1) + H2(o["conn_Sz"], o["conn_Sp"], Sz1, Sp1),
        "conn_Sz": kron(identity(my_block), Sz1),
        "conn_Sp": kron(identity(my_block), Sp1),
    }
    
    return Enlarged_Block(length=block.length + 1,  basis_size=block.basis_size * model1_d, operator_dict=enlarge_block_dict)


def rotate_and_truncate(operator, transformation_matrix):
    """
    Transforms the operator to the new (possibly truncated) basis given by transformation matrix.
    :param operator: operator to transform (np.array).
    :param transformation_matrix:   transformation matrix (np.array).
    :return: transformed matrix.
    """
    return transformation_matrix.conjugate().transpose().dot(operator.dot(transformation_matrix))


def single_dmrg_step(sys, env, m):
    """
    Single step of DMRG.
    :param sys: system (dictionary)
    :param env: environment (dictionary)
    :param m: max number of states (int)
    :return: new block and energy (tuple)
    """
    assert is_valid_block(sys)
    assert is_valid_block(env)

    sys_enl = enlarge_block(sys)
    if sys is env:
        env_enl = sys_enl
    else:
        env_enl = enlarge_block(env)
    
    assert is_valid_enlarged_block(sys_enl)
    assert is_valid_enlarged_block(env_enl)

    m_sys_enl = sys_enl.basis_size
    m_env_enl = env_enl.basis_size
    sys_enl_op = sys_enl.operator_dict
    env_enl_op = env_enl.operator_dict
    super_block_hamiltonian = kron(sys_enl_op["H"], identity(m_env_enl)) + kron(identity(m_sys_enl), env_enl_op["H"]) \
                              + H2(sys_enl_op["conn_Sz"], sys_enl_op["conn_Sp"], env_enl_op["conn_Sz"], env_enl_op["conn_Sp"])

    (energy,), psi0 = eigsh(super_block_hamiltonian, k=1, which="SA")

    psi0 = psi0.reshape([sys_enl.basis_size, -1], order="C")
    rho = np.dot(psi0, psi0.conjugate().transpose())

    e_values, e_vectors = np.linalg.eigh(rho)
    newmman = sum([-(x * np.log2(x)) for x in e_values])
    possible_eigenstates = []
    for e_value, e_vector in zip(e_values, e_vectors.transpose()):
        possible_eigenstates.append((e_value, e_vector))
    possible_eigenstates.sort(reverse=True, key=lambda x: x[0])

    my_m = min(len(possible_eigenstates), m)
    transformation_matrix = np.zeros((sys_enl.basis_size, my_m), dtype="d", order='F')
    for i, (e_value, e_vector) in enumerate(possible_eigenstates[:my_m]):
        transformation_matrix[:, i] = e_vector

    truncation_error = 1 - sum(x[0] for x in possible_eigenstates[:my_m])

    new_operator_dict = {}
    for name, op in sys_enl.operator_dict.items():
        new_operator_dict[name] = rotate_and_truncate(op, transformation_matrix)
    
    new_block = Block(length=sys_enl.length,
                      basis_size=my_m,
                      operator_dict=new_operator_dict)
    return new_block, energy, truncation_error, newmman


def infinite_system_algorithm(L, m):
    """
    Implementation of simple infinite system dmrg.
    :param L: long of system (int)
    :param m: max number of states (int)
    :return: None.
    """
    block = initial_block
    folder_1 = open("datam="+str(m)+".txt", "a")
    folder_2 = open("longm="+str(m)+".txt", "a")

    while 2 * block.length < L:
        folder_2.write(str(block.length * 2 + 2) + "\n")
        block, energy, truncation_error, newmman = single_dmrg_step(block, block, m=m)
        folder_1.write(str(energy) +"\t"+str(newmman)+"\t" + str(truncation_error) + "\n")
        print(str(energy) +"\t"+str(newmman)+"\t" + str(truncation_error) + "\n")
    folder_1.close()
    folder_2.close()


def main():
    np.set_printoptions(precision=10, suppress=True, threshold=10000, linewidth=300)
    for k in [12]:
        infinite_system_algorithm(L=500, m=k)


if __name__ == '__main__':
    main()

