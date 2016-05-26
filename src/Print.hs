{- Print
by Gregory Schwartz

-- | Collects the functions pertaining to the printing of the BLOSUM matrix
-}

{-# LANGUAGE OverloadedStrings #-}

module Print
    ( printBlosum
    , printBlosumCSV
    ) where

-- Standard
import Data.Maybe
import qualified Data.Map.Strict as Map

-- Cabal
import qualified Data.Text as T
import TextShow

-- Local
import Types

-- | Print the BLOSUM matrix as a dataframe
printBlosum :: Blosum -> T.Text
printBlosum = T.append header . T.append "\n" . body . unBlosum
  where
    header = "source,destination,value"
    body   = T.unlines
           . concatMap Map.elems
           . Map.elems
           . Map.mapWithKey
             ( \(AA x) -> Map.mapWithKey
                          (\ (AA y) (BlosumVal v)
                          -> T.intercalate "," [ T.singleton x
                                               , T.singleton y
                                               , showt v
                                               ]
                          )
             )

-- | Print the BLOSUM matrix as a csv matrix according to a certain order
printBlosumCSV :: [AA] -> Blosum -> T.Text
printBlosumCSV order = T.append header . T.append "\n" . body . unBlosum
  where
    header     = T.intercalate "," . ("" :) . map (T.singleton . unAA) $ order
    body m     = T.unlines
               . map ( \x -> T.append (flip T.snoc ',' . T.singleton . unAA $ x)
                           . fromMaybe ""
                           . flip Map.lookup (printMap m)
                           $ x
                     )
               $ order
    printMap :: Map.Map AA (Map.Map AA BlosumVal) -> Map.Map AA T.Text
    printMap   = Map.map ( \v -> T.intercalate ","
                               . map (flip lookZero v)
                               $ order
                         )
    lookZero k = fromMaybe "" . fmap (showt . unBlosumVal) . Map.lookup k
