{- Print
by Gregory Schwartz

-- | Collects the functions pertaining to the printing of the BLOSSUM matrix
-}

{-# LANGUAGE OverloadedStrings #-}

module Print
    ( printBlossum
    , printBlossumCSV
    ) where

-- Standard
import Data.Maybe
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

-- | Print the BLOSSUM matrix as a csv matrix according to a certain order
printBlossumCSV :: [AA] -> Blossum -> T.Text
printBlossumCSV order = T.append header . T.append "\n" . body . unBlossum
  where
    header     = T.intercalate "," . ("" :) . map (T.singleton . unAA) $ order
    body m     = T.unlines
               . map ( \x -> T.append (flip T.snoc ',' . T.singleton . unAA $ x)
                           . fromMaybe ""
                           . flip Map.lookup (printMap m)
                           $ x
                     )
               $ order
    printMap :: Map.Map AA (Map.Map AA BlossumVal) -> Map.Map AA T.Text
    printMap   = Map.map ( \v -> T.intercalate ","
                               . map (flip lookZero v)
                               $ order
                         )
    lookZero k = fromMaybe "" . fmap (showt . unBlossumVal) . Map.lookup k
