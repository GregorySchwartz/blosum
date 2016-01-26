{- Print
by Gregory Schwartz

-- | Collects the functions pertaining to the printing of the BLOSSUM matrix
-}

{-# LANGUAGE OverloadedStrings #-}

module Print
    ( printBlossum
    ) where

-- Standard
import qualified Data.Map.Strict as Map

-- Cabal
import qualified Data.Text as T
import TextShow

-- Local
import Types

-- | Print the BLOSSUM matrix as a dataframe
printBlossum :: Blossum -> T.Text
printBlossum = T.append header . T.append "\n" . body . unBlossum
  where
    header = "source,destination,value"
    body   = T.unlines
           . concatMap Map.elems
           . Map.elems
           . Map.mapWithKey
             ( \(AA x) -> Map.mapWithKey
                          (\ (AA y) (BlossumVal v)
                          -> T.intercalate "," [ T.singleton x
                                               , T.singleton y
                                               , showt v
                                               ]
                          )
             )
