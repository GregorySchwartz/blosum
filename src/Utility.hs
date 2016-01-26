{- Utility
by Gregory Schwartz

-- | Collects the functions pertaining to generic helpful functions
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE OverloadedStrings #-}

module Utility
    ( getSeq
    , sortTuple
    , nub'
    , groupBlocks
    , cutEnds
    ) where

-- Standard
import Data.List
import qualified Data.Set as Set
import qualified Data.Sequence as Seq

-- Cabal
import qualified Data.Text as T
import Data.Fasta.Text

-- Local
import Types

-- | Faster nub
nub' :: (Ord a) => [a] -> [a]
nub' = Set.toList . Set.fromList

-- | Extract sequence from a fasta sequence
getSeq :: FastaSequence -> [AA]
getSeq = map AA . T.unpack . fastaSeq

-- | Sort a tuple
sortTuple :: (Ord a) => (a, a) -> (a, a)
sortTuple t@(!x, !y)
    | x <= y    = t
    | otherwise = (y, x)

-- | Group together fasta entires by a field
groupBlocks :: Maybe Field -> [FastaSequence] -> [[FastaSequence]]
groupBlocks (Just field) = groupBy
                           (\x y -> getField field x == getField field y)
groupBlocks Nothing      = (:[])

-- | Get a field from a fasta header
getField :: Field -> FastaSequence -> T.Text
getField (Field x) = (!! x) . T.splitOn "|" . fastaHeader

-- | Cut off the ends of a sequence
cutEnds :: Seq.Seq a -> Seq.Seq a
cutEnds = Seq.reverse . Seq.drop 1 . Seq.reverse . Seq.drop 1
