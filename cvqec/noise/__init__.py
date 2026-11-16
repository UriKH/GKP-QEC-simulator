"""Noise channels. Add one by subclassing Channel and decorating it with @register_channel."""

from cvqec.noise.base import CHANNELS, Channel, Compose, Identity, make_channel, register_channel
from cvqec.noise.channels import Dephasing, GaussianDisplacement, PhotonLoss

__all__ = ["CHANNELS", "Channel", "Compose", "Dephasing", "GaussianDisplacement", "Identity", "PhotonLoss",
           "make_channel", "register_channel"]
