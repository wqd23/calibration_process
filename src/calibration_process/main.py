from lib_reader import single_read03b


def main():
    sci, tel = single_read03b(
        "data/03B/raw_data/20210501_tempbias_03B/-10C_265_4m_rundata2021-05-01-14-27-37.dat",
        "data/03B/raw_data/20210501_tempbias_03B/-10C_265_4m_scienceConfig2021-05-01-14-27-37.json",
    )
    print(sci.keys())


if __name__ == "__main__":
    main()
