    def inverse_region_RNA(self, interval, tss_file):
        st, en = self.get_interval(interval, return_ids=True)
        old_length = len(self.chr_data[interval.chr])
        tss_file = pd.read_csv(tss_file, encoding='utf-8', sep='\t')
        # Добавить колонку с направлением транскрипции путем слияния данных RNAseq и tss файла:
        self.chr_data[interval.chr] = pd.merge(self.chr_data[interval.chr],
                                               tss_file.loc[:, ['Strand', 'Gene stable ID']], how="left",
                                               left_on="gene", right_on="Gene stable ID")
        # Поиск генов, затронутых инверсией. Создание списка их индексов:
        drop_indices = list(np.where(
            (((self.chr_data[interval.chr].Strand == 1) & ((self.chr_data[interval.chr].start < interval.start) &
                                                           ((self.chr_data[interval.chr].start + 2000 * self.chr_data[
                                                               interval.chr].Strand) > interval.start)) |
              ((self.chr_data[interval.chr].start < interval.end) &
               ((self.chr_data[interval.chr].start + 2000 * self.chr_data[interval.chr].Strand) > interval.end)))) |
            (((self.chr_data[interval.chr].Strand == -1) & ((self.chr_data[interval.chr].end > interval.start) &
                                                            ((self.chr_data[interval.chr].end + 2000 * self.chr_data[
                                                                interval.chr].Strand) < interval.start)) |
              ((self.chr_data[interval.chr].end > interval.end) &
               ((self.chr_data[interval.chr].end + 2000 * self.chr_data[interval.chr].Strand) < interval.end)))))[0])

        debug = len(self.chr_data[interval.chr])
        # print(drop_indices, '  dropped')
        # print('i.end', interval.end)
        # y = (self.chr_data[interval.chr].iloc[10:20, [self.chr_data[interval.chr].columns.get_loc("chr"),
        #                                               self.chr_data[interval.chr].columns.get_loc("start"),
        #                                               self.chr_data[interval.chr].columns.get_loc("end"),
        #                                               self.chr_data[interval.chr].columns.get_loc("sigVal"),
        #                                               self.chr_data[interval.chr].columns.get_loc("Strand"),
        #                                               self.chr_data[interval.chr].columns.get_loc("gene")]])
        # print(y)

        # inv_indices = list(self.chr_data[interval.chr].index[(self.chr_data[interval.chr].start > interval.start) &
        #                                                      (self.chr_data[interval.chr].end < interval.end)])
        # if len(inv_indices) > 0: print("inv_indices 1\n", inv_indices)
        # inv_data = self.chr_data[interval.chr].loc[inv_indices]  # inversed genes as
        # if len(inv_data) > 0: print("inv data 1\n", inv_data)

        # # exit()
        # # for i, x in self.chr_data[interval.chr].iloc[st:en].iterrows():
        # #     print(x["end"])
        # # exit()
        # # for x in self.chr_data[interval.chr].iloc[st:en + 1, [self.chr_data[interval.chr].columns.get_loc("end")]]
        # #    print("x,type(x)\n",x,type(x))
        # # Change coordinates of inversed genes
        # for i, x in self.chr_data[interval.chr].iloc[st:en].iterrows():
        #     if x["end"] > (interval.start + interval.len / 2):
        #         inv_data["start"] = interval.start + (interval.end - x["end"])
        #     else:
        #         inv_data["start"] = interval.end - (x["end"] - interval.start)
        # print("inv_data 2\n", inv_data)
        # exit()
        # inv_data["end"] = [interval.start + (interval.end - x) if x > (interval.start + interval.len / 2)
        #                    else interval.end - (x - interval.start) for x in self.chr_data[interval.chr].iloc[st:en + 1,
        #                                                                      self.chr_data[
        #                                                                          interval.chr].columns.get_loc(
        #                                                                          "start")]]
        # print("dup data 3\n", inv_data)
        # exit()
        # print(interval.start + interval.len / 2,'interval.start + interval.len / 2')
        """новые координаты для инвертируемых генов"""
        starts = [interval.start + (interval.end - x) if x > (interval.start + interval.len / 2)
                  else interval.end - (x - interval.start) for x in
                  self.chr_data[interval.chr].iloc[st:en + 1, self.chr_data[interval.chr].columns.get_loc("end")]]
        # print(starts, 'starts')
        ends = [interval.start + (interval.end - x) if x > (interval.start + interval.len / 2)
                else interval.end - (x - interval.start) for x in
                self.chr_data[interval.chr].iloc[st:en + 1, self.chr_data[interval.chr].columns.get_loc("start")]]
        # print(ends, 'ends')
        # exit()
        if len(starts) > 0:
            self.chr_data[interval.chr].iloc[st:en + 1, self.chr_data[interval.chr].columns.get_loc("start")] = starts
        if len(ends) > 0:
            self.chr_data[interval.chr].iloc[st:en + 1, self.chr_data[interval.chr].columns.get_loc("end")] = ends

        self.chr_data[interval.chr].sort_values(by="start", inplace=True)

        # Delete genes affected by invertion:
        if len(drop_indices) > 0:
            self.chr_data[interval.chr].drop(drop_indices, inplace=True)

        # Set new indices according new "start" and "end":
        self.chr_data[interval.chr].set_index(self.chr_data[interval.chr].apply(
            lambda x: pd.Interval(x.start, x.end, closed="both"), axis="columns"), inplace=True)
        self.chr_data[interval.chr].drop(['Strand', 'Gene stable ID'], axis=1, inplace=True)
        assert old_length - debug == 0
        # print(old_length, len(self.chr_data[interval.chr]))

        # # print('self.chr_data[interval.chr]','\n' ,self.chr_data[interval.chr])
        # # Set new indices according new "start" and "end":
        # self.chr_data[interval.chr].set_index(self.chr_data[interval.chr].apply(
        #     lambda x: pd.Interval(x.start, x.end, closed="both"), axis="columns"), inplace=True)
        # self.chr_data[interval.chr].drop(['Strand', 'Gene stable ID'], axis=1, inplace=True)
        # # assert len(self.chr_data[interval.chr]) - debug == old_length


"""начало шевеления клеткой мозга 24.03.22"""

""""необходимое место для полимеразы = 2000 от старта транскрипции гена в сторону экспресии. Если попадёт дальше, то считаем, что 3D организация не меняется
нужно посмотреть, каккие гены пересекаются с границами перестройки
потом поменять координаты из 
starts and ends from chipseq