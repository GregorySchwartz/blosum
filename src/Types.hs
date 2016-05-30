{- Types
by Gregory Schwartz

-- | Collects the functions pertaining to the types used in the program
-}

{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Types where

-- Standard
import qualified Data.Map.Strict as Map
import qualified Data.Sequence as Seq

-- Cabal
import Data.Fasta.Text (FastaSequence)

newtype AA        = AA { unAA :: Char } deriving (Eq, Ord, Show, Read)
newtype Nuc       = Nuc Char deriving (Eq, Ord, Show, Read)
newtype Field     = Field Int deriving (Eq)
newtype Frequency = Frequency Double
                    deriving (Eq, Ord, Num, Enum, Show, Read, Fractional)
newtype Position  = Position Int deriving (Eq, Ord, Num, Enum, Show, Read)
newtype Identity  = Identity Double deriving (Eq, Ord, Num, Show, Read)
newtype BlosumVal = BlosumVal { unBlosumVal :: Int }
                    deriving (Eq, Ord, Num, Enum, Show, Read)

newtype ClusterMap          =
    ClusterMap { unClusterMap :: Map.Map Int (Seq.Seq FastaSequence) }
    deriving (Eq, Ord, Show)
newtype ClusterFrequencyMap =
    ClusterFrequencyMap
    { unClusterFrequencyMap :: (Map.Map Position (Seq.Seq (AA, Frequency)))
    }
    deriving (Eq, Ord, Show)
newtype AAMap               = AAMap
                              { unAAMap :: Map.Map AA (Map.Map AA Frequency) }
                              deriving (Eq, Ord, Show)
newtype BlockMap            = BlockMap { unBlockMap :: AAMap }
                              deriving (Eq, Ord, Show)
newtype FrequencyMap        = FrequencyMap { unFrequencyMap :: AAMap }
                              deriving (Eq, Ord, Show)
newtype Blosum              =
    Blosum { unBlosum :: (Map.Map AA (Map.Map AA BlosumVal)) }
    deriving (Eq, Ord, Show)

instance Monoid AAMap where
    mempty                      = AAMap Map.empty
    mappend (AAMap x) (AAMap y) = AAMap $ Map.unionWith (Map.unionWith (+)) x y
