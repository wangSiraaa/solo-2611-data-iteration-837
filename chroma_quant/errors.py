"""统一异常类型。"""


class ChromaError(Exception):
    """本包所有异常的基类。"""


class ValidationError(ChromaError):
    """输入数据不满足约束（类型、单调性、非负性等）。"""


class LimitError(ChromaError):
    """超过题面给定的规模上限。"""


class DomainError(ChromaError):
    """查询/积分落在有效域之外：本服务拒绝外推。"""


class NotFoundError(ChromaError):
    """run、版本、产品、方法或分析不存在。"""
