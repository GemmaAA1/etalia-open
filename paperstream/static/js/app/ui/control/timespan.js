define(['jquery', 'app/util/utils'], function ($, Utils) {

    var Timespan = function (options) {
        this.config = $.extend({
            debug: false,
            element: '#timespan',
            toggle: '#toggle-timespan',
            selection: '#timespan-selection'
        }, options);

        this.$element = $(this.config.element);
        this.$toggle = $(this.config.toggle);
        this.$selection = $(this.config.selection);
    };

    Timespan.prototype.init = function() {
        var that = this,
            $body = $('body');

        this.$toggle.on('click', function(e) {
            Utils.toggleClass(that.$element, 'opened');

            e.preventDefault();
            return false;
        });
        // Close on click out
        $(window).on('click', function(e) {
            if (0 == $(e.target).closest(that.config.element).length) {
                that.$element.removeClass('opened');
            }
        });

        this.$element.on('click', '.choices a', function(e) {
            that.setValue($(e.target).closest('a').data('timespan'));

            that.$toggle.trigger('click');

            $body.trigger('etalia.control.timespan.change');

            e.preventDefault();
            return false;
        });

        return this;
    };

    Timespan.prototype.setValue = function(value) {
        var text = (function(s) {
            switch (s) {
                case 7 :  return 'W';
                case 30 : return '1m';
                case 60 : return '2m';
            }
            throw 'Unexpected timespan value';
        })(parseInt(value));

        this.$selection.html(text).data('value', value);

        return this;
    };

    Timespan.prototype.getValue = function() {
        return this.$selection.data('value');
    };

    return Timespan;
});
