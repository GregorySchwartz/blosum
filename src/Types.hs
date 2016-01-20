{- Types
by Gregory Schwartz

-- | Collects the functions pertaining to the types used in the program
-}

{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Types where

newtype Identity = Identity Double deriving (Eq, Ord, Num, Show, Read)
