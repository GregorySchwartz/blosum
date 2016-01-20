{- Cluster
by Gregory Schwartz

-- | Collects the functions pertaining to the clustering of sequences by hamming
distance identity
-}

module Cluster
    ( clusterIdentity
    ) where

-- Standard
import qualified Data.Set.Monad as Set

-- Cabal
import Data.Fasta.Text
import Math.Diversity.Diversity (hamming)

-- | Takes in an identity and fasta sequences and returns the sequences
-- grouped together by hamming distance identity
clusterIdentity :: Identity
                -> [FastaSequence]
                -> Map.Map Int (Set.Set FastaSequence)
clusterIdentity identity = 

-- | Keep comparing clusters until no more fusions (no change in size) make
-- sense
clusterIdentityGo :: Bool
                  -> Identity
                  -> [Set.Set FastaSequence]
                  -> [Set.Set FastaSequence]
compareGroups identity xs
    | 
  where
    before = Map.size xs

-- | Either the groups are joined by identity or not
compareGroups :: Identity
              -> Seq.Set FastaSequence
              -> Seq.Set FastaSequence
              -> Maybe Seq.Set FastaSequence
compareGroups identity xs ys
    | any $ (\x y -> getIdentity x y >= identity) <$> xs <*> ys =
        xs Set.union ys
    | otherwise                                                 = Nothing

-- | Get the identity between two sequences
getIdentity :: FastaSequence -> FastaSequence -> Identity
getIdentity xs ys = (* 100)
                  . (/ T.length xs)
                  . hamming (fastaSeq xs)
                  $ (fastaSeq ys)
