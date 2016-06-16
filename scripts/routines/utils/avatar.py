# -*- coding: utf-8 -*-
from __future__ import unicode_literals, absolute_import
import os
from random import randint, seed
from io import BytesIO
from PIL import Image, ImageDraw, ImageFont
import tempfile


class AvatarGenerator(object):
    FONT_COLOR = (255, 255, 255)
    MIN_RENDER_SIZE = 512

    def generate(self, size, string, filetype="JPEG"):
        """
            Generates a squared avatar with random background color.
            :param size: size of the avatar, in pixels
            :param string: string to be used to print text and seed the random
            :param filetype: the file format of the image (i.e. JPEG, PNG)
        """
        render_size = max(size, self.MIN_RENDER_SIZE)
        image = Image.new('RGB', (render_size, render_size),
                          self._background_color(string))
        draw = ImageDraw.Draw(image)
        font = self._font(render_size)
        text = self._text(string)
        draw.text(self._text_position(render_size, text, font),
                  text,
                  fill=self.FONT_COLOR,
                  font=font)
        # stream = BytesIO()
        f = tempfile.NamedTemporaryFile()
        image = image.resize((size, size), Image.ANTIALIAS)
        filename = '{0}.{1}'.format(f.name, filetype)
        image.save(filename, format=filetype, optimize=True)
        return filename

    @staticmethod
    def _background_color(s):
        """
            Generate a random background color.
            Brighter colors are dropped, because the text is white.
            :param s: Seed used by the random generator
            (same seed will produce the same color).
        """
        seed(s)
        r = v = b = 255
        while r + v + b > 255*2:
            r = randint(0, 255)
            v = randint(0, 255)
            b = randint(0, 255)
        return (r, v, b)

    @staticmethod
    def _font(size):
        """
            Returns a PIL ImageFont instance.
            :param size: size of the avatar, in pixels
        """
        path = os.path.join(os.path.dirname(__file__), '../data',
                            "Inconsolata.otf")
        return ImageFont.truetype(path, size=int(0.8 * size))

    @staticmethod
    def _text(string):
        """
            Returns the text to draw.
        """
        if len(string) == 0:
            return "#"
        else:
            return string[0].upper()

    @staticmethod
    def _text_position(size, text, font):
        """
            Returns the left-top point where the text should be positioned.
        """
        width, height = font.getsize(text)
        left = (size - width) / 2.0
        # I just don't know why 5.5, but it seems to be the good ratio
        top = (size - height) / 5.5
        return left, top