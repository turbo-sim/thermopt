"""Verify actual optional media encoding and the manager's public frame rate."""

from pathlib import Path

import imageio.v2 as imageio
import numpy as np
from PIL import Image
import pytest

import thermopt as th


@pytest.fixture
def animation_manager(tmp_path):
    manager = th.ThermodynamicCycleOptimization.__new__(
        th.ThermodynamicCycleOptimization
    )
    manager.out_dir = str(tmp_path)
    directory = tmp_path / "optimization"
    directory.mkdir()
    for index in range(3):
        pixels = np.zeros((64, 96, 3), dtype=np.uint8)
        pixels[:, :, index] = 200
        imageio.imwrite(directory / f"iteration_{index:03d}.png", pixels)
    return manager


@pytest.mark.parametrize("format", ["gif", "mp4"])
def test_animation_preserves_frame_count_order_and_requested_frame_rate(
    animation_manager, format
):
    manager = animation_manager
    manager.create_animation(format=format, fps=4)
    filename = Path(manager.out_dir) / f"optimization_animation.{format}"
    if format == "gif":
        with Image.open(filename) as image:
            assert image.n_frames == 3
            for index in range(image.n_frames):
                image.seek(index)
                assert image.info["duration"] == 250
                assert np.argmax(np.asarray(image.convert("RGB"))[10, 10]) == index
    else:
        with imageio.get_reader(filename) as reader:
            assert reader.get_meta_data()["fps"] == pytest.approx(4)
            assert reader.count_frames() == 3
            for index in range(3):
                assert np.argmax(reader.get_data(index)[10, 10]) == index


@pytest.mark.parametrize(
    "format,fps", [("png", 1), ("gif", 0), ("mp4", -1), ("both", float("nan"))]
)
def test_animation_rejects_unsupported_requests_before_writing(
    animation_manager, format, fps
):
    with pytest.raises(ValueError, match="format|fps"):
        animation_manager.create_animation(format=format, fps=fps)
    assert not list(Path(animation_manager.out_dir).glob("optimization_animation.*"))
