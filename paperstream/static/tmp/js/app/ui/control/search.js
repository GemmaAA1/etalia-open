define(['jquery', 'app/util/utils'], function ($, Utils) {

    var Search = function (options) {
        this.config = $.extend({
            debug: false,
            element: '#search',
            toggle: '#toggle-search',
            close: '#close-search',
            input: '#search-input'
        }, options);

        this.$element = $(this.config.element);

        this.$toggle = this.$element.find(this.config.toggle);
        this.$close  = this.$element.find(this.config.close);
        this.$input  = this.$element.find(this.config.input);

        this.keyUpTimeout = null;
    };

    Search.prototype.init = function() {
        var that = this,
            $body = $('body');

        this.$toggle.on('click', function(e) {
            Utils.toggleClass($search, 'opened');

            e.preventDefault();
            return false;
        });

        this.$close.on('click', function(e) {
            Utils.toggleClass($search, 'opened');

            e.preventDefault();
            return false;
        });

        this.$input
            .on('focus', function(e) {
                e.stopPropagation();
                $(e.delegateTarget).parents('form').addClass('active');
            })
            .on('blur', function(e) {
                e.stopPropagation();
                $(e.delegateTarget).parents('form').removeClass('active');
            })
            .on('keyup', function(e) {
                if (that.keyUpTimeout) {
                    clearTimeout(that.keyUpTimeout);
                }
                var code = e.keyCode || e.which;
                if (code == 13) { // Enter pressed
                    $body.trigger('etalia.control.search.change');
                } else {
                    that.keyUpTimeout = setTimeout(function() {
                        $body.trigger('etalia.control.search.change');
                    }, 1000);
                }
            });

        return this;
    };

    Search.prototype.setValue = function(value) {
        this.$input.val(value);

        return this;
    };

    Search.prototype.getValue = function() {
        return this.$input.val();
    };

    return Search;
});
