import shutil
import requests


def download_file(url, write_path):
    with open(write_path, "wb") as f:
        response = requests.get(url, stream=True)
        total_length = response.headers.get("content-length")

        assert total_length is not None

        progress = 0
        total_length = int(total_length)
        for data in response.iter_content(131072):
            progress += len(data)
            f.write(data)

            bar = int(80 * progress / total_length)
            print(
                f"\rDownloading: {progress/1024**2:.2f}/{total_length/1024**2:.2f}MB [{'='*bar}{' '*(80-bar)}]",
                end="",
            )
    print("")


def unpack_file(path_to_file, write_path):
    shutil.unpack_archive(path_to_file, write_path)


if __name__ == "__main__":
    # DOI: https://zenodo.org/doi/10.5281/zenodo.11118021
    out_file = "tmp/scs_analysis_data.zip"

    print("Downloading Dataset (https://zenodo.org/doi/10.5281/zenodo.11118021)")

    download_file(
        "https://zenodo.org/records/11118022/files/scs_analysis_data.zip?download=1",
        out_file,
    )
    print("Download Complete!")

    print("Unpacking...")
    unpack_file(out_file, "data/")
    print("Finished!")
