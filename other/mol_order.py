# coding=utf-8

from enum import Enum


class MolOrder(Enum):
    """Configure how to resort molecules

    KEEP_ORDER : Keep molecules order as is.
    MOL_ORDER : Keep molecule types order as is. Molecules with same type will
                grouped together. Molecules in each type groups will keep
                their original order.
    MOL_ORDER_ALPHA : Molecules with same type will be grouped together, while
                      molecule types will be reordered in alphabetic order.
                      Molecules in each type groups will keep their original order.
    """
    KEEP_ORDER = 0
    MOL_ORDER = 1
    MOL_ORDER_ALPHA = 2
