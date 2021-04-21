module Dupuy_Galichon_Matching
export hello, domath

"""
    hello(who::String)

Return "Hello, `who`".
"""
hello(who::String) = "Hello, $who"

"""
    domath(x::Number, y::Number)

Return `x + y`.
"""
domath(x::Number, y::Number) = x + y

end
